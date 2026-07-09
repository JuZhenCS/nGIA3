/*
clustering.hpp
输入文件数据内容:
  uint32_t  序列数
  vector<size_t>  packed数据偏移 readsCount
  vector<size_t>  fasta数据偏移  readsCount
  分组数据  (128 * readsCount条记录)
    vector<uint32_t> 序列对应的代表序列 readsCount
    ... (共128组)
    vector<uint32_t> 序列对应的代表序列 readsCount
  packed数据 (readsCount条记录)
    uint32_t 数据长度 1
    uint32_t 压缩数据 (数据长度+31)/32*5
  fasta数据 (readsCount条记录)
    string 序列名 1
    string 序列 1
长度<0xFFFF 数量<0x7FFFFFFF
用法:
ngia2 clustering -p packed文件 -r result文件 -i identity值
2025-08-04 by 鞠震
*/

#ifndef CLUSTERING_HPP
#define CLUSTERING_HPP

#include <iostream>  // cout
#include <vector>  // vector
#include <cstdint>  // uint32_t
#include <fstream>  // ifstream
#include <cstring>  // memset
#include <cuda.h>  // cuda
#include "timer.hpp"  // timer
#include "parser.hpp"  // parser

namespace Clustering {  // 命名空间

//--------数据--------//
struct Option {  // 参数
  std::string packedFile;  // packed文件
  std::string resultFile;  // 结果文件
  uint32_t identity;  // 相似度
};

struct Data{  // 数据
  uint32_t seqCount;  // 序列
  size_t *packedOffsets;  // packed数据偏移
  std::vector<size_t> fastaOffsets;  // fasta数据偏移
  std::vector<uint32_t> groups;  // 分组信息
  uint32_t *packedReads;  // packed数据
};

//--------函数--------//
// init 初始化 ok
void init(int argc, char **argv, Option &option) {
  {  // 解析命令行
    Parser::Parser parser;
    parser.add(true, "-p", "packed file (string)", "");
    parser.add(true, "-r", "result file (string)", "");
    parser.add(true, "-i", "identity (1-99)", "");
    parser.parse(argc, argv);
    option.packedFile = parser.getValue<std::string>("-p");
    option.resultFile = parser.getValue<std::string>("-r");
    option.identity = parser.getValue<uint32_t>("-i");
  }
  bool checkPassed = true;  // 每一项校验都通过才进行下一步
  {  // 校验参数
    std::ifstream packedFile(option.packedFile);
    if (!packedFile.good()) {  // 没有输入文件
      std::cout << option.packedFile << " does not exists.\n";
      checkPassed = false;
    }
    if (option.identity<1 || option.identity>99) {  // 相似度溢出
      std::cout << "identity should be 1-99\n";
      checkPassed = false;
    }
  }
  cudaDeviceProp prop;  // 显卡属性
  {  // 配置显卡 export CUDA_VISIBLE_DEVICES=0 指定GPU1
    cudaSetDevice(0);
    cudaGetDeviceProperties(&prop, 0);
    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);  // 共享内存 缓存优先
    cudaDeviceSynchronize();
    if (cudaGetLastError() != cudaSuccess) {  // GPU状态校验
      std::cout << cudaGetErrorString(cudaGetLastError()) << "\n";
      checkPassed = false;
    }
  }
  if (!checkPassed) { exit(0); }
  std::cout << "use gpu:\t" << prop.name << "\n";
  std::cout << "packed file:\t" << option.packedFile << "\n";
  std::cout << "result file:\t" << option.resultFile << "\n";
  std::cout << "identity:\t" << option.identity << "\n";
}

// readData 读取数据
void readData(const Option &option, Data &inputData) {
  uint32_t seqCount = 0;  // 序列数
  std::ifstream packedFile(option.packedFile);
  packedFile.read((char*)&seqCount, sizeof(uint32_t));  // 序列数
  inputData.seqCount = seqCount;
  cudaMallocManaged(&inputData.packedOffsets, sizeof(size_t)*seqCount);
  packedFile.read((char*)inputData.packedOffsets, sizeof(size_t)*seqCount);
  inputData.fastaOffsets.resize(seqCount, 0);  // fasta数据偏移
  packedFile.read((char*)inputData.fastaOffsets.data(),sizeof(size_t)*seqCount);
  inputData.groups.resize(128*seqCount, 0);  // 分组信息
  packedFile.read((char*)inputData.groups.data(),sizeof(uint32_t)*128*seqCount);
  size_t tempLength = 0;  // 打包数据长度
  tempLength = inputData.fastaOffsets[0] - inputData.packedOffsets[0];
  cudaMallocManaged(&inputData.packedReads, tempLength);  // 打包数据
  packedFile.read((char*)inputData.packedReads, tempLength);  // 读打包数据
  tempLength = inputData.packedOffsets[0];  // 把文件偏移转为内存数据偏移
  for (uint32_t i=0; i<seqCount; i++) {  // 字节位置转为uint32_t偏移
    inputData.packedOffsets[i] -= tempLength;
    inputData.packedOffsets[i] /= sizeof(uint32_t);  // char转uint32_t
  }
  std::cout << "reads count:\t" << seqCount << "\n";  // 序列数
  std::cout << "longest:\t" << inputData.packedReads[0] << "\n";  // 最长
  tempLength = inputData.packedOffsets[seqCount-1];  // 最短序列的位置
  std::cout << "shortest:\t" << inputData.packedReads[tempLength] << "\n";
}

// 不要改内外循环 寄存器使用会变少
// 不要数据预取 或操作指针 用线程数掩盖延迟
// kernel_dynamic 动态规划
__global__ void __launch_bounds__(64, 1)  // maxThreadblock, minBlock/SM
kernel_dynamic(const uint32_t *reads, const size_t *offsets,
const uint32_t *jobs, const uint32_t jobCount, uint32_t *clusters,
const float threshold) {
  uint32_t index = blockDim.x*blockIdx.x+threadIdx.x; // 线程编号
  if (index >= jobCount) { return; }  // 超出范围
  const uint32_t* __restrict__ represent = reads+offsets[jobs[index*2]];  //代表
  const uint32_t* __restrict__ target = reads+offsets[jobs[index*2+1]];  // 目标
  const uint32_t lengthRep = represent[0];  // 代表序列长度
  const uint32_t lengthTar = target[0];  // 目标序列长度
// 修改
  if (lengthTar < lengthRep*threshold) { return; }  // 代表短于目标
// 修改
  uint32_t lines[2048];  // 每行结果 别全赋初值 也别动态申请 开销太大
  memset(lines, 0xFF, (lengthRep+31)/32*sizeof(uint32_t));  // 0是匹配 1是不匹配
  uint32_t Rows[5] = {0};  // 从代表取的32个碱基/氨基酸 行
  uint32_t Cols[5] = {0};  // 从目标取的32个碱基/氨基酸 列
  uint32_t matchs[29] = {0};  // 匹配的碱基/氨基酸 1是匹配 0是不匹配
  // 左偏移要对齐32 否则+32以后，会超出长度范围，导致漏比对
  const uint32_t lsft = ((uint32_t)ceil(lengthTar-lengthTar*threshold)+31)&~31;
  // 右偏移要加32 因为是分块计算，一块长是32，按右下角算右偏移
  const uint32_t rsft = (uint32_t)ceil(lengthRep-lengthTar*threshold)+32;
  for (uint32_t i=0; i<lengthTar; i+=32) {  // 遍历列
    uint32_t carrys = 0;  // 进位
    #pragma unroll  // 展开循环得到列
    for (uint32_t e=0; e<5; e++) { Cols[e] = target[1+i/32*5+e]; }  // 取列
    const uint32_t jstart = max((int32_t)i-(int32_t)lsft, (int32_t)0);  // 开始
    const uint32_t jend = min((int32_t)i+(int32_t)rsft, (int32_t)lengthRep);
    for (uint32_t j=jstart; j<jend; j+=32) {  // 遍历行
      #pragma unroll  // 展开循环得到行
      for(uint32_t e=0; e<5; e++) { Rows[e] = represent[1+j/32*5+e]; }
      for (uint32_t k=1; k<29; k++) {  // 预生成match 0跟谁都不匹配
        uint32_t match = 0xFFFFFFFF;
        #pragma unroll  // 展开循环得到单个match
        for (uint32_t e=0; e<5; e++) { match &= Rows[e]^0xFFFFFFFF+(k>>e&1); }
        matchs[k] = match;
      }
      uint32_t row = lines[j/32];  // 上一行结果
      for (uint32_t k=0; k<32; k++) {  // 32*32的核心
        uint32_t order = 0;
        for (uint32_t e=0; e<5; e++) { order += (Cols[e]>>k&1)<<e; }
        uint32_t match = matchs[order];  // 匹配上的碱基/氨基酸
        uint32_t carry = carrys & 1;  // 进位
        uint32_t term0 = row & match;
        uint32_t term1 = row & (~match);
        uint32_t carryRow = row + carry;
        carry = carryRow < row;  // 是否发生进位
        carryRow += term0;
        carry |= carryRow < term0;  // 是否发生进位
        row = carryRow | term1;
        carrys = (carrys>>1) + (carry<<31);  // 写回进位
      }
      lines[j/32] = row;  // 写回结果
    }
  }
  // 统计结果
  uint32_t lcs = 0;  // 最长公共子序列长度
  for (uint32_t i=0; i<lengthRep; i+=32) { lcs += 32-__popc(lines[i/32]); }
// 修改
  if (lcs > (lengthRep+lengthTar-lcs)*threshold) {  // 不用优化 没第二个分支了
  // if (lcs > lengthRep*threshold) {  // 不用优化 没第二个分支了
  // if (lcs > lengthTar*threshold) {  // 不用优化 没第二个分支了
// 修改
    clusters[jobs[index*2+1]] = jobs[index*2];  // 写入聚类结果
  }
}

// clustering 利用局部敏感哈希快速聚类
void clustering(const Option &option, const Data &inputData,
std::vector<uint32_t> &results) {
  // 中间数据
  uint32_t seqCount = inputData.seqCount;  // 序列数
  uint32_t *clusters = NULL;  // 聚类结果 clusters[i]是i的聚类结果
  cudaMallocManaged(&clusters, sizeof(uint32_t)*seqCount);
  memset(clusters, 0xFF, sizeof(uint32_t)*seqCount);  // 最大值是未聚类
  uint32_t *jobs = NULL;  // 比对任务 代表序列+目标序列 2个uint32_t一组
  cudaMallocManaged(&jobs, sizeof(uint32_t)*seqCount*2);
  memset(jobs, 0, sizeof(uint32_t)*seqCount*2);  // 赋初始方便debug
  const float threshold = option.identity/100.0f;  // 序列相似度


  // // 聚类过程
  for (uint32_t loop=0; loop<128; loop++) {  // 重复分组比对过程
    std::cout << "\rclustering:\t" << loop+1 << "/128" << std::flush;
    uint32_t jobCount = 0;
    for (uint32_t i=0; i<seqCount; i++) {  // 生成任务
      uint32_t represent = inputData.groups[loop*seqCount+i];  // 代表序列
      uint32_t target = i;  // 目标序列
      if (clusters[target]!=0xFFFFFFFF || target==represent) { continue; }
      jobs[jobCount*2] = represent;  // 代表序列
      jobs[jobCount*2+1] = target;  // 目标序列
      jobCount += 1;  // 任务数+1
    }
    kernel_dynamic<<<(jobCount+63)/64, 64>>>(
      inputData.packedReads, inputData.packedOffsets,
      jobs, jobCount, clusters, threshold
    );
    cudaDeviceSynchronize();  // 等数据传输完成
  }


  // // 聚类过程 两两比对1
  // for (uint32_t loop=0; loop<seqCount; loop++) {  // 重复分组比对过程
  //   std::cout << "\rclustering:\t" << loop+1 << "/" << seqCount << std::flush;
  //   uint32_t jobCount = 0;
  //   for (uint32_t i=loop+1; i<seqCount; i++) {  // 生成任务
  //     uint32_t represent = loop;  // 代表序列
  //     uint32_t target = i;  // 目标序列
  //     if (clusters[i]==0xFFFFFFFF) {  // 没有聚类
  //       jobs[jobCount*2] = represent;  // 代表序列
  //       jobs[jobCount*2+1] = target;  // 目标序列
  //       jobCount += 1;  // 任务数+1
  //     }
  //   }
  //   kernel_dynamic<<<(jobCount+63)/64, 64>>>(
  //     inputData.packedReads, inputData.packedOffsets,
  //     jobs, jobCount, clusters, threshold
  //   );
  //   cudaDeviceSynchronize();  // 等数据传输完成
  // }


  // // 聚类过程 两两比对2,  跳过代表序列 是原始的贪婪增量方法
  // for (uint32_t loop=0; loop<seqCount; loop++) {  // 重复分组比对过程
  //   std::cout << "\rclustering:\t" << loop+1 << "/" << seqCount << std::flush;
  //   if (clusters[loop] != 0xFFFFFFFF) { continue; }  // 跳过非代表序列
  //   uint32_t jobCount = 0;
  //   for (uint32_t i=loop+1; i<seqCount; i++) {  // 生成任务
  //     uint32_t represent = loop;  // 代表序列
  //     uint32_t target = i;  // 目标序列
  //     if (clusters[i]==0xFFFFFFFF) {  // 没有聚类
  //       jobs[jobCount*2] = represent;  // 代表序列
  //       jobs[jobCount*2+1] = target;  // 目标序列
  //       jobCount += 1;  // 任务数+1
  //     }
  //   }
  //   kernel_dynamic<<<(jobCount+63)/64, 64>>>(
  //     inputData.packedReads, inputData.packedOffsets,
  //     jobs, jobCount, clusters, threshold
  //   );
  //   cudaDeviceSynchronize();  // 等数据传输完成
  // }


  std::cout << std::endl;
  // 收尾
  results.resize(seqCount, 0);  // 聚类结果
  for (uint32_t i=0; i<seqCount; i++) { results[i] = clusters[i]; }  // 复制结果
  cudaFree(inputData.packedOffsets);  // 释放内存
  cudaFree(inputData.packedReads);  // 释放内存
  cudaFree(clusters);  // 释放内存
  cudaFree(jobs);  // 释放内存
}

// saveResult 保存聚类结果
void saveResult(const Option &option, const Data &inputData,
  const std::vector<uint32_t> &results) {
  const uint32_t seqCount = results.size();  // 序列数
  // 整理聚类结果
  std::vector<std::vector<uint32_t>> finalGroups(0);  // 最终分组
  finalGroups.resize(seqCount, std::vector<uint32_t>(0));  // 初始化
  for (uint32_t i=0; i<seqCount; i++) {  // 遍历所有序列
    if (results[i]==0xFFFFFFFF) {  // 是代表序列
      finalGroups[i].push_back(i);  // 加入自己
    } else {  // 不是代表序列
      uint32_t represent = results[i];  // 直到找到根代表
      while (results[represent]!=0xFFFFFFFF) { represent = results[represent]; }
      finalGroups[represent].push_back(i);  // 加入代表序列的组
    }
  }
  // 写入文件
  std::ifstream fastaFile(option.packedFile, std::ios::in);  // 输入文件
  std::ofstream resultFile(option.resultFile, std::ios::out);  // 输出文件
  std::string line = "";  // 临时字符串
  uint32_t groupCount = 0;
  for (uint32_t i=0; i<seqCount; i++) {  // 遍历所有序列
    if (finalGroups[i].size() == 0) { continue; }  // 跳过空组
    uint32_t represent = finalGroups[i][0];  // 代表序列
    fastaFile.seekg(inputData.fastaOffsets[represent], std::ios::beg);  // 定位
    std::getline(fastaFile, line);  // 读取标题
    resultFile << line << std::endl;  // 写入标题
    std::getline(fastaFile, line);  // 读取序列
    resultFile << line << std::endl;  // 写入序列
    for (uint32_t j=1; j<finalGroups[i].size(); j++) {  // 遍历所有成员
      uint32_t member = finalGroups[i][j];  // 成员序列
      fastaFile.seekg(inputData.fastaOffsets[member], std::ios::beg);  // 定位
      std::getline(fastaFile, line);  // 读取标题
      resultFile << "  "+line << std::endl;  // 写入标题
    }
    groupCount += 1;  // 分组数+1
  }
  fastaFile.close();  // 关闭输入
  resultFile.close();  // 关闭输出
  std::cout << "group count:\t" << groupCount << std::endl;
}

//--------主函数--------//
void clustering(int argc, char **argv) {
  Timer::Timer timer;
  timer.printStamp();

  Option option = { packedFile:"", resultFile:"", identity:0 };
  init(argc, argv, option);  // 解析命令行
  Data inputData = {  // 输入数据
    seqCount:0, packedOffsets:NULL, fastaOffsets:std::vector<size_t>(0),
    groups:std::vector<uint32_t>(0), packedReads:NULL
  };
  readData(option, inputData);  // 读取数据
  std::vector<uint32_t> results(0);  // 聚类结果
  clustering(option, inputData, results);  // 聚类
  saveResult(option, inputData, results);  // 保存结果

  timer.printStamp();
  timer.printDuration();
}
}  // namespace Clustering

#endif  // CLUSTERING_HPP

// 优化
// 比对算法优化
// cudaMemAdvise(reads, position-offsets[0], cudaMemAdviseSetReadMostly, 0);
// 数据预取
// -maxrregcount 56 --resource-usage
// cudaMemPrefetchAsync(remains, sizeof(uint32_t)*remainCount, cudaCpuDeviceId,
// 0); cudaMemPrefetchAsync(remains, sizeof(uint32_t)*remainCount, 0);
// cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);  // 共享内存变缓存
// __restrict__
// 常量内存
// 显存内计算remains
// 少量序列用batch

// 4090的SM参数:
// 1536个线程
// 48个warp
// 24个block
// 64K个寄存器
// 八个warp就能隐藏延迟了，四发射，64线程足够
