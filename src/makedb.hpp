/*
makeDB.hpp
输出文件数据内容:
  uint32_t         序列数         1
  vector<size_t>   packed数据偏移  seqCount
  vector<size_t>   fasta数据偏移   seqCount
  分组数据 (128*readsCount条记录)
    vector<uint32_t> 序列对应的代表序列 readsCount
    ... (共128组)
    vector<uint32_t> 序列对应的代表序列 readsCount
  packed数据 (readsCount条记录)
    uint32_t 数据长度 1
    uint32_t 压缩数据 (数据长度+31)/32*5
  fasta数据 (readsCount条记录)
    string 序列名 1
    string 序列 1
长度<=0xFFFF 数量<=0x7FFFFFFF
用法:
ngia2 makeDB -f fasta文件 -p packed文件
2026-01-15 by 鞠震
*/

// 序列长 < 0xFFFF     比对时line  < 2048
// 序列数 < 0x7FFFFFFF 比对时线程数 < int范围

#ifndef MAKEDB_HPP
#define MAKEDB_HPP

#include <iostream>  // std::cout
#include <fstream>  // fstream
#include <unordered_map>  // unordered_map
#include <algorithm>  // stable_sort
#include <omp.h>  // openmp
#include "timer.hpp"  // timer
#include "parser.hpp"  // parser
#include <random>  // mt19937_64


namespace MakeDB {  // 命名空间

//--------数据--------//
struct Option {  // 输入选项
  std::string fastaFile;  // fasta文件
  std::string packedFile;  // packed文件
};
struct Sequence {  // 记录序列
  std::string name;  // 序列名字
  std::string read;  // 序列数据
};

// 修改

// const std::unordered_map<char, uint32_t> transTable = {  // 转码表
//   {'a',  1}, {'b',  2}, {'c',  3}, {'d',  4}, {'e',  5}, {'f',  6}, {'g',  7},
//   {'A',  1}, {'B',  2}, {'C',  3}, {'D',  4}, {'E',  5}, {'F',  6}, {'G',  7},
//   {'h',  8}, {'i',  9}, {'j', 10}, {'k', 11}, {'l', 12}, {'m', 13}, {'n', 14},
//   {'H',  8}, {'I',  9}, {'J', 10}, {'K', 11}, {'L', 12}, {'M', 13}, {'N', 14},
//   {'o', 15}, {'p', 16}, {'q', 17}, {'r', 18}, {'s', 19}, {'t', 20}, {'u', 21},
//   {'O', 15}, {'P', 16}, {'Q', 17}, {'R', 18}, {'S', 19}, {'T', 20}, {'U', 21},
//   {'v', 22}, {'w', 23}, {'x', 24}, {'y', 25}, {'z', 26}, {'*', 27},
//   {'V', 22}, {'W', 23}, {'X', 24}, {'Y', 25}, {'Z', 26}, {'-', 28}
// };

// const std::unordered_map<char, uint32_t> transTable = {  // linclust的转码表
//   {'a',  1}, {'b',  2}, {'c',  3}, {'d',  4}, {'e',  5}, {'f',  6}, {'g',  7},
//   {'A',  1}, {'B',  2}, {'C',  3}, {'D',  4}, {'E',  5}, {'F',  6}, {'G',  7},
//   {'h',  8}, {'i',  9}, {'j', 10}, {'k', 11}, {'l', 12}, {'m', 12}, {'n',  4},
//   {'H',  8}, {'I',  9}, {'J', 10}, {'K', 11}, {'L', 12}, {'M', 12}, {'N',  4},
//   {'o', 13}, {'p', 14}, {'q',  5}, {'r', 11}, {'s',  1}, {'t', 17}, {'u', 18},
//   {'O', 13}, {'P', 14}, {'Q',  5}, {'R', 11}, {'S',  1}, {'T', 17}, {'U', 18},
//   {'v',  9}, {'w', 19}, {'x', 20}, {'y',  6}, {'z', 21},
//   {'V',  9}, {'W', 19}, {'X', 20}, {'Y',  6}, {'Z', 21}
// };

// const std::unordered_map<char, uint32_t> transTable = {  // 精简转码表
//   {'a',  1}, {'c',  2}, {'d',  3}, {'e',  4}, {'f',  5}, {'g',  6}, {'h',  7},
//   {'A',  1}, {'C',  2}, {'D',  3}, {'E',  4}, {'F',  5}, {'G',  6}, {'H',  7},
//   {'i',  8}, {'k',  9}, {'l', 10}, {'m', 10}, {'n',  3}, {'o', 13}, {'p', 14},
//   {'I',  8}, {'K',  9}, {'L', 10}, {'M', 10}, {'N',  3}, {'O', 13}, {'P', 14},
//   {'q',  4}, {'r',  9}, {'s',  1}, {'t', 15}, {'u', 19}, {'v',  8}, {'w', 21},
//   {'Q',  4}, {'R',  9}, {'S',  1}, {'T', 14}, {'U', 19}, {'V',  8}, {'W', 21},
//   {'y',  5}, {'Y',  5},
// };

const std::unordered_map<char, uint32_t> transTable = {  // 精简转码表
  {'a',  1}, {'c',  2}, {'d',  3}, {'e',  4}, {'f',  5}, {'g',  6}, {'h',  7},
  {'A',  1}, {'C',  2}, {'D',  3}, {'E',  4}, {'F',  5}, {'G',  6}, {'H',  7},
  {'i',  8}, {'k',  9}, {'l', 10}, {'m', 10}, {'n',  3}, {'o', 11}, {'p', 12},
  {'I',  8}, {'K',  9}, {'L', 10}, {'M', 10}, {'N',  3}, {'O', 11}, {'P', 12},
  {'q',  4}, {'r',  9}, {'s',  1}, {'t', 13}, {'u', 14}, {'v',  8}, {'w', 15},
  {'Q',  4}, {'R',  9}, {'S',  1}, {'T', 13}, {'U', 14}, {'V',  8}, {'W', 15},
  {'y',  5}, {'Y',  5}
};

// const std::unordered_map<char, uint32_t> transTable = {  // BLOSUM 转码表 2-3分
//   {'a',  1}, {'b',  2}, {'c',  3}, {'d',  4}, {'e',  4}, {'f',  6}, {'g',  7},
//   {'A',  1}, {'B',  2}, {'C',  3}, {'D',  4}, {'E',  4}, {'F',  6}, {'G',  7},
//   {'h',  8}, {'i',  9}, {'j', 10}, {'k', 11}, {'l',  9}, {'m',  9}, {'n', 14},
//   {'H',  8}, {'I',  9}, {'J', 10}, {'K', 11}, {'L',  9}, {'M',  9}, {'N', 14},
//   {'o', 15}, {'p', 16}, {'q',  4}, {'r', 11}, {'s', 19}, {'t', 20}, {'u', 21},
//   {'O', 15}, {'P', 16}, {'Q',  4}, {'R', 11}, {'S', 19}, {'T', 20}, {'U', 21},
//   {'v',  9}, {'w',  6}, {'x', 24}, {'y',  6}, {'z', 26}, {'*', 27},
//   {'V',  9}, {'W',  6}, {'X', 24}, {'Y',  6}, {'Z', 26}, {'-', 28}
// };

// const std::unordered_map<char, uint32_t> transTable = {  // BLOSUM 转码表 3分
//   {'a',  1}, {'b',  2}, {'c',  3}, {'d',  4}, {'e',  5}, {'f',  6}, {'g',  7},
//   {'A',  1}, {'B',  2}, {'C',  3}, {'D',  4}, {'E',  5}, {'F',  6}, {'G',  7},
//   {'h',  8}, {'i',  9}, {'j', 10}, {'k', 11}, {'l', 12}, {'m', 13}, {'n', 14},
//   {'H',  8}, {'I',  9}, {'J', 10}, {'K', 11}, {'L', 12}, {'M', 13}, {'N', 14},
//   {'o', 15}, {'p', 16}, {'q', 17}, {'r', 18}, {'s', 19}, {'t', 20}, {'u', 21},
//   {'O', 15}, {'P', 16}, {'Q', 17}, {'R', 18}, {'S', 19}, {'T', 20}, {'U', 21},
//   {'v',  9}, {'w', 23}, {'x', 24}, {'y',  6}, {'z', 26}, {'*', 27},
//   {'V',  9}, {'W', 23}, {'X', 24}, {'Y',  6}, {'Z', 26}, {'-', 28}
// };

// 修改

//--------函数--------//
// init 初始化
void init(int argc, char **argv, Option &option) {
  {  // 解析命令行
    Parser::Parser parser;  // 解析器
    parser.add(true, "-f", "fasta file (string)", "");  // fasta
    parser.add(true, "-p", "packed file (string)", "");  // packed
    if (!parser.parse(argc, argv)) { exit(0); }  // 解析命令行
    option.fastaFile = parser.getValue<std::string>("-f");  // fasta文件
    option.packedFile = parser.getValue<std::string>("-p");  // packed文件
  }
  bool checkPassed = true;  // 校验是否通过
  {  // 校验参数
    std::ifstream fastaFile(option.fastaFile);  // fasta文件
    if (!fastaFile.good()) {  // 文件不存在
      std::cout << option.fastaFile << " does not exist.\n";
      checkPassed = false;
    }
    std::string line = "";  // 一行数据
    getline(fastaFile, line);
    if (line.back() == '\r') {  // 行尾有\r
      std::cout << "Please remove the \\r at the end of lines.\n";
      checkPassed = false;
    }
  }
  if (!checkPassed) { exit(0); }  // 校验
  std::cout << "fasta:\t" << option.fastaFile << "\n";
  std::cout << "packed:\t" << option.packedFile << "\n";
  std::cout << "thread:\t" << omp_get_max_threads() << "\n";
}

// readData 读数据
void readData(const std::string &fastaFile, std::vector<Sequence> &sequences) {
  Sequence sequence = {name:"", read:""};  // 一条序列
  std::string line = "";  // 一行数据
  uint32_t seqCount = 0;  // 序列数
  std::ifstream fileIn(fastaFile);  // fasta文件
  std::cout << "read:\t" << std::flush;  // 进度条
  while (fileIn.peek()!=EOF && sequences.size()<=0x7FFFFFFF) {  // 读文件
    sequence.name = "";  // 序列名字
    sequence.read = "";  // 序列数据
    getline(fileIn, sequence.name);  // 读序列名
    while (fileIn.peek()!=EOF && fileIn.peek()!='>') {  // 序列数<=2147483647
      getline(fileIn, line);
      sequence.read += line;
    }
    if (sequence.read.size() <= 0xFFFF) { sequences.push_back(sequence); }
    if (seqCount%(1024*1024) == 0) std::cout << "." << std::flush;
    seqCount+=1;
  }
  fileIn.close();  // 关闭文件
  sequences.shrink_to_fit();  // 省点内存
  std::stable_sort(sequences.begin(), sequences.end(),  // 排序
    [](const Sequence &a, const Sequence &b) {
      return a.read.size() > b.read.size(); });
  std::cout << " finish\n";
  std::cout << "count:\t" << seqCount << "\n";
}

// 打包数据
void packData(const std::vector<Sequence> &sequences,
std::vector<std::vector<uint32_t>> &packedReads) {
  uint32_t seqCount = sequences.size();  // 序列数
  packedReads.resize(seqCount, std::vector<uint32_t>(0));  // 打包后数据
  std::cout << "pack:\t" << std::flush;  // 打印进度
  #pragma omp parallel for schedule(dynamic, 16) // 并行打包
  for (uint32_t i=0; i<seqCount; i++) {  // 遍历序列
    packedReads[i].assign(1+(sequences[i].read.size()+31)/32*5, 0);  // 打包数据
    uint32_t length = 0;  // 已压缩长度
    for (char base : sequences[i].read) {  // 遍历碱基
      std::unordered_map<char, uint32_t>::const_iterator it;  // 查找转码表
      it = transTable.find(base);
      uint32_t pack = it==transTable.end() ? 0 : it->second;  // 转码
      for (uint32_t e=0; e<5; e++) {  // 打包
        packedReads[i][1+length/32*5+e] += (pack>>e&1)<<(length%32);
      }
      length += 1;
    }
    packedReads[i][0] = length;  // 长度
    if (i%(1024*1024) == 0) std::cout << "." << std::flush;  // 打印进度
  }
  std::cout << " finish\n";  // 打印进度
}

// 签名哈希
inline uint32_t bioWangHash(uint32_t x, uint32_t seed) {
  x = x ^ seed;
  x = (x ^ 61) ^ (x >> 16);
  x = x + (x << 3);
  x = x ^ (x >> 4);
  x = x * 0x27d4eb2d;
  x = x ^ (x >> 15);
  // 额外的雪崩效应
  x = (x ^ (x >> 8)) * 0x9e3779b1;
  x = x ^ (x >> 14);
  return x;
}

uint32_t lcg_hash(uint32_t x) {
    // 一个常用的 LCG 参数，来自 Numerical Recipes
    x = x * 1664525u + 1013904223u;
    return x;
}

uint32_t multiplicative_hash(uint32_t x) {
    // 0x9e3779b1 是常用的黄金比例常数
    return x * 0x9e3779b1u;
}

uint32_t xorshift_hash(uint32_t x) {
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return x;
}

// uint32_t get_seed_for_round(int round) {
//     // 使用一个固定种子初始化生成器，确保每次运行结果可复现
//     static std::mt19937_64 rng(0x9e3779b97f4a7c15ULL); 
//     static std::uniform_int_distribution<uint32_t> dist(0, UINT32_MAX);
    
//     // 为每一轮生成一个不同的随机种子
//     return dist(rng); 
// }

uint32_t get_seed_for_round(int round) {
    uint64_t base_seed = 0x9e3779b97f4a7c15ULL; // 一个固定的黄金比例常数
    // 将轮次索引与固定基数混合，并截断为32位
    return (uint32_t)((base_seed ^ (round * 0x9e3779b9)) & 0xFFFFFFFF);
}


// 序列签名
void signSequence(const std::vector<std::vector<uint32_t>> &packedReads,
std::vector<uint32_t> &signsTable) {
  const uint32_t seqCount = packedReads.size();  // 序列数
  signsTable.assign(seqCount*128, 0xFFFFFFFF);  // 先填充
  std::cout << "sign:\t." << std::flush;
  #pragma omp parallel for schedule(dynamic, 16)  // 并行签名
  for (uint32_t i=0; i<seqCount; i++) {
    const std::vector<uint32_t> packed = packedReads[i];  // 打包数据
    const uint32_t length = packed[0];  // 序列长度
// 修改
    uint32_t kmer = 0;  // kmer
// 修改
    for (uint32_t j=0; j<length; j++) {  // 遍历kmer
      kmer <<= 4;  // 腾位置
      for (uint32_t e=0; e<4; e++) {kmer+=((packed[1+j/32*5+e]>>(j%32))&1)<<e;}
// 修改
      kmer &= (1<<25)-1;  // 截取固定位
      if (j<6 && j<(length-1)) {continue;}  // 前后5个碱基不签名
// 修改
      for (uint32_t k=0; k<128; k++) {  // 签名128次
        uint32_t seed = k * 0x9e3779b1 + 0x85ebca6b;  // 随机种子
        uint32_t sign = bioWangHash(kmer, seed);  // 签名哈希
        // uint32_t sign = xorshift_hash(kmer+seed);  // 签名哈希
        if (sign < signsTable[i*128+k]) { signsTable[i*128+k] = sign; }  // 更新
      }
    }
  }
  std::cout << " finish\n";  // 打印进度
}

// // 序列签名，只用稀有kmer签名
// void signSequence(const std::vector<std::vector<uint32_t>> &packedReads,
// std::vector<uint32_t> &signsTable) {
//   const uint32_t seqCount = packedReads.size();  // 序列数
//   signsTable.assign(seqCount*128, 0xFFFFFFFF);  // 先填充
//   std::cout << "sign:\t." << std::flush;
//   #pragma omp parallel for schedule(dynamic, 16)  // 并行签名
//   for (uint32_t i=0; i<seqCount; i++) {
//     // 找到序列的稀有k-mer
//     std::unordered_map <uint32_t, uint32_t> kmerCount;  // kmer计数
//     const std::vector<uint32_t> packed = packedReads[i];  // 打包数据
//     const uint32_t length = packed[0];  // 序列长度
//     uint32_t kmer = 0;  // kmer
//     if (length <= 5) {  // 长度不足一个签名，有多少签多少
//       for (uint32_t j=0; j<length; j++) {  // 遍历kmer
//         kmer <<= 5;  // 腾位置
//         for (uint32_t e=0; e<5; e++) {kmer+=((packed[1+j/32*5+e]>>(j%32))&1)<<e;}
//       }
//       kmerCount[kmer] = 1;  // 记录到kmerCounnt
//     } else {  // 长度足够一个签名，需要截取固定长度
//       for (uint32_t j=0; j<4; j++) {  // kmer基础
//         kmer <<= 5;  // 腾位置
//         for (uint32_t e=0; e<5; e++) {kmer+=((packed[1+j/32*5+e]>>(j%32))&1)<<e;}
//       }
//       for (uint32_t j=4; j<length; j++) {  // 遍历kmer
//         kmer <<= 5;  // 腾位置
//         for (uint32_t e=0; e<5; e++) {kmer+=((packed[1+j/32*5+e]>>(j%32))&1)<<e;}
//         kmer &= (1<<30)-1;  // 截取固定位
//         if (kmerCount.find(kmer) != kmerCount.end()) {  // kmerCount中找到
//           kmerCount[kmer] += 1;  // 计数加1
//         } else {  // 没找到
//           kmerCount[kmer] = 1;  // 记录到kmerCounnt
//         }
//       }
//     }
//     for (auto it = kmerCount.begin(); it != kmerCount.end(); it++) {  // 遍历kmerCount
//       if (it->second > 1) {  // 稀有kmer
//         uint32_t kmer = it->first;  // 一个签名
//         for (uint32_t k=0; k<128; k++) {  // 签名128次
//           uint32_t seed = k * 0x9e3779b1 + 0x85ebca6b;  // 随机种子
//           uint32_t sign = bioWangHash(kmer, seed);  // 签名哈希
//           signsTable[i*128+k] = std::min(signsTable[i*128+k], sign);  // 签名
//         }
//       }
//     }
//   }
//   std::cout << " finish\n";  // 打印进度
// }


// 序列预分组
// void groupSequence(const std::vector<uint32_t> &signsTable,
// std::vector<uint32_t> &groupTable) {
//   std::cout << "group:\t." << std::flush;
//   const uint32_t seqCount = signsTable.size()/128;  // 序列数
//   std::vector<std::vector<std::pair<int32_t, uint32_t>>> signs(0);  // 转置签名
//   signs.assign(128,std::vector<std::pair<int32_t, uint32_t>>(seqCount,{0,0}));
//   #pragma omp parallel for schedule(dynamic, 1)  // 并行转置
//   for (uint32_t i=0; i<128; i++) {  // 遍历序列
//     for (uint32_t j=0; j<seqCount; j++) {  // 遍历签名
//       signs[i][j].first = signsTable[j*128+i];  // 转置签名
//       signs[i][j].second = j;  // 序列编号
//     }
//     std::stable_sort(signs[i].begin(), signs[i].end(),
//       [](const std::pair<int32_t, uint32_t> &a,
//       const std::pair<int32_t, uint32_t> &b) { return a.first < b.first; }
//     );
//   }
//   groupTable.assign(128*seqCount, 0);  // 分组表
//   #pragma omp parallel for schedule(dynamic, 1)  // 并行回填
//   for (uint32_t i=0; i<128; i++) {  // 把分组结果写回
//     uint32_t sign = signs[i][0].first;  // 签名
//     uint32_t represent = signs[i][0].second;  // 代表序列
//     for (uint32_t j=0; j<seqCount; j++) {  // 遍历序列
//       if (signs[i][j].first != sign) {  // 签名不同
//         sign = signs[i][j].first;  // 更新签名
//         represent = signs[i][j].second;  // 更新代表序列
//       }
//       groupTable[i*seqCount+signs[i][j].second] = represent;  // 写入分组表
//     }
//   }
//   std::cout << " finish\n";  // 打印进度
// }

// 修改
// 序列预分组
void groupSequence(const std::vector<uint32_t> &signsTable,
std::vector<uint32_t> &groupTable) {
  std::cout << "group:\t." << std::flush;
  const uint32_t seqCount = signsTable.size()/128;  // 序列数
  std::vector<std::vector<std::pair<int32_t, uint32_t>>> signs(0);  // 转置签名
  signs.assign(128,std::vector<std::pair<int32_t, uint32_t>>(seqCount,{0,0}));
  #pragma omp parallel for schedule(dynamic, 1)  // 并行转置
  for (uint32_t i=0; i<128; i++) {  // 遍历序列
    for (uint32_t j=0; j<seqCount; j++) {  // 遍历签名
      signs[i][j].first = signsTable[j*128+i];  // 转置签名
      signs[i][j].second = j;  // 序列编号
    }
    std::stable_sort(signs[i].begin(), signs[i].end(),
      [](const std::pair<int32_t, uint32_t> &a,
      const std::pair<int32_t, uint32_t> &b) { return a.first < b.first; }
    );
  }
  groupTable.assign(128*seqCount, 0);  // 分组表
  #pragma omp parallel for schedule(dynamic, 1)  // 并行回填
  for (uint32_t i=0; i<128; i++) {  // 把分组结果写回
    uint32_t sign = signs[i][0].first+1;  // 签名要与第一个不同
    uint32_t represent = 0;  // 代表序列
    for (uint32_t j=0; j<seqCount; j++) {  // 遍历序列
      if (signs[i][j].first != sign) {  // 签名不同
        sign = signs[i][j].first;  // 更新签名
        represent = signs[i][j].second;  // 更新代表序列
      } else {  // 签名相同
        represent = signs[i][j-1].second;  // 更新代表序列
      }
      groupTable[i*seqCount+signs[i][j].second] = represent;  // 写入分组表
    }
  }
  std::cout << " finish\n";  // 打印进度
}
// 修改

// 写入数据
void writeFile(const std::string &packedFile,
const std::vector<Sequence> &sequences,
const std::vector<std::vector<uint32_t>> &packedReads,
const std::vector<uint32_t> groupTable) {
  std::cout << "save:\t." << std::flush;
  const uint32_t seqCount = sequences.size();  // 序列数
  std::vector<size_t> packedOffsets(seqCount, 0);  // packed偏移
  std::vector<size_t> fastaOffsets(seqCount, 0);  // fasta偏移
  size_t offset = sizeof(uint32_t)+sizeof(size_t)*seqCount*2;  // 临时偏移
  offset += sizeof(uint32_t)*128*seqCount;  // 还有groups
  for (uint32_t i=0; i<seqCount; i++) {  // 计算packed偏移
    packedOffsets[i] = offset;
    offset += sizeof(uint32_t)*(packedReads[i].size());
  }
  for (uint32_t i=0; i<seqCount; i++) {  // fasta偏移
    fastaOffsets[i] = offset;
    offset += sequences[i].name.size()+sequences[i].read.size()+2;  // 包括换行
  }
  std::ofstream fileOut(packedFile);  // 输出文件 写入结果
  fileOut.write((char*)&seqCount, sizeof(uint32_t));  // 写入序列数
  fileOut.write((char*)packedOffsets.data(), sizeof(size_t)*seqCount);  // 偏移
  fileOut.write((char*)fastaOffsets.data(), sizeof(size_t)*seqCount);  // 偏移
  fileOut.write((char*)groupTable.data(), sizeof(uint32_t)*groupTable.size());
  for (uint32_t i=0; i<seqCount; i++) {  // 写入packed数据
    uint32_t length = packedReads[i].size();  // 长度
    fileOut.write((char*)packedReads[i].data(), sizeof(uint32_t)*length);
  }
  for (uint32_t i=0; i<seqCount; i++) {  // 写入fasta数据
    std::string line = sequences[i].name+"\n"+sequences[i].read+"\n";
    fileOut.write(line.c_str(), line.size());
  }
  fileOut.close();  // 关闭文件
  std::cout << " finish\n";  // 打印进度
}

//--------主函数--------//
void makeDB(int argc, char **argv) {
  Timer::Timer timer;  // 计时器
  timer.printStamp();  // 时间戳

  Option option = {fastaFile:"", packedFile:""};  // 选项
  init(argc, argv, option);  // 初始化
  std::vector<Sequence> sequences(0);  // 序列集合
  readData(option.fastaFile, sequences);  // 读入数据
  std::vector<std::vector<uint32_t>> packedReads(0);  // 打包后数据
  packData(sequences, packedReads);  // 打包数据
  std::vector<uint32_t> signsTable(0);  // 签名表
  signSequence(packedReads, signsTable);  // 序列签名
  std::vector<uint32_t> groupTable(0);  // 分组表
  groupSequence(signsTable, groupTable);  // 序列分组
  writeFile(option.packedFile, sequences, packedReads, groupTable);  // 写入文件

  timer.printStamp();  // 时间戳
  timer.printDuration();  // 耗时
}

}  // namespace MakeDB
#endif  // MAKEDB_HPP

// 需要测试以下内容对精度和速度的影响：
// 1. kmer长度
// 2. 哈希后收束到序列数范围  nmi大幅降低了 速度也变慢了
// 3. r和b的值 r是1已经是精度最高的方法了
// 4. 满32个再生成kmer，还是5个就生成kmer；5个就要生成，否则就糊了
// 5. 切kmer的时候，切的大一点？
// 6. 试试加权minhash
// 7. 分组的时候，考虑序列长度
