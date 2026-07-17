/*
makeDB.hpp
输出文件数据内容:
  uint32_t         序列数           1
  uint32_t[seqCount]  nameLengths  序列名字长度
  uint32_t[seqCount]  readLengths  序列数据长度
  vector<size_t>   packed数据偏移  seqCount
  vector<size_t>   fasta数据偏移   seqCount
  签名数据 (128*readsCount条记录)
    vector<uint32_t> 序列对应的签名 readsCount
    ... (共128组)
  packed数据 (readsCount条记录)
    uint32_t 数据长度 1
    uint32_t 压缩数据 (数据长度+31)/32*4
  fasta数据 (readsCount条记录)
    string 序列名含换行
    string 序列含换行
长度<=0xFFFF 数量<=0x7FFFFFFF
用法:
ngia2 makeDB -f fasta文件 -p packed文件
2026-07-09 by 鞠震
*/

// 序列长 < 0xFFFF     比对时line  < 2048
// 序列数 < 0x7FFFFFFF 比对时线程数 < int范围

#ifndef MAKEDB_HPP
#define MAKEDB_HPP

#include <iostream>  // cout
#include <fstream>   // ifstream/fstream
#include <vector>    // vector
#include <string>    // string
#include <algorithm> // stable_sort
#include <array>     // array
#include <cstdint>   // uint32_t
#include <omp.h>     // openmp
#include "timer.hpp" // timer
#include "parser.hpp" // parser

namespace MakeDB {  // 命名空间

//--------------------------------  数据结构  --------------------------------//
struct Option {  // 输入选项
  std::string fastaFile;  // fasta文件
  std::string packedFile; // packed文件
};

struct seqIndex {  // 序列索引
  size_t fastaOffset;  // fasta文件偏移
  size_t packedOffsetPack;  // packed文件中 打包偏移
  size_t packedOffsetSeq;  // packed文件中 序列偏移
  uint32_t nameLength;  // 序列名字长度
  uint32_t readLength;  // 序列数据长度
  uint32_t seqLength;  // 序列记录长度 名字+数据+换行
};

//--------------------------------  功能函数  --------------------------------//
// init初始化
Option init(int argc, char** argv) {
  Option option = {fastaFile:"", packedFile:""};  // 初始化为空
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
    // 1. 文件存在
    std::ifstream fastaFile(option.fastaFile);  // fasta文件
    if (!fastaFile.good()) {  // 文件不存在
      std::cout << option.fastaFile << " does not exist.\n";
      checkPassed = false;
    }
    // 2. 行尾没有\r，只\n换行
    std::string line;  // 一行数据
    getline(fastaFile, line);
    if (!line.empty() && line.back() == '\r') {  // 行尾有\r
      std::cout << "Please remove the \\r at the end of lines.\n";
      checkPassed = false;
    }
    // 3. 文件以换行结尾
    fastaFile.seekg(0, std::ios::end);
    if (fastaFile.tellg() > 0) {
      fastaFile.seekg(-1, std::ios::end);
      char lastChar;
      fastaFile.read(&lastChar, 1);
      if (lastChar != '\n') {
        std::cout << "File must end with \\n.\n";
        checkPassed = false;
      }
    }
  }
  if (!checkPassed) { exit(0); }  // 校验
  std::cout << "fasta:\t" << option.fastaFile << "\n";  // 输入数据
  std::cout << "packed:\t" << option.packedFile << "\n";  // 打包输出
  std::cout << "thread:\t" << omp_get_max_threads() << "\n";  // 可用线程数
  return option;
}

// makeIndex 第一遍读文件，构建索引
std::vector<seqIndex> makeIndex(const std::string& fastaFile) {
  std::vector<seqIndex> indices;  // 序列索引
  uint32_t seqCount = 0;  // 扫描的序列数（含超长丢弃的）
  {  // 1. 从文件读序列索引
    std::cout << "read:\t." << std::flush;  // 进度条
    std::ifstream fileIn(fastaFile);  // fasta文件
    std::string line = "";  // 一行数据
    size_t fastaOffset = 0;  // fasta偏移
    uint32_t nameLength = 0;  // 序列名字长度
    uint32_t readLength = 0;  // 序列数据长度
    uint32_t seqLength = 0;  // 序列长度
    while (fileIn.peek() != EOF && seqCount < 0x7FFFFFFF) {
      fastaOffset = fileIn.tellg();  // fasta偏移
      getline(fileIn, line);  // 读 >name 行
      nameLength = line.size();  // 序列名字长度
      seqLength = line.size() + 1;  // 包括换行符
      readLength = 0;  // 序列数据长度
      while (fileIn.peek() != EOF && fileIn.peek() != '>') {
        getline(fileIn, line);  // 序列数据
        readLength += line.size();  // 序列数据长度
        seqLength += line.size() + 1;  // 包括换行符
      }
      seqCount += 1;  // 文件中的序列总数
      if (readLength < 0xFFFF) {  // 只保留长度小于65536的序列
        indices.push_back({
          fastaOffset: fastaOffset,
          packedOffsetPack: 0,
          packedOffsetSeq: 0,
          nameLength: nameLength,
          readLength: readLength,
          seqLength: seqLength
        });
      }
      if (seqCount % (1024 * 1024) == 0) { std::cout << "." << std::flush; }
    }
    std::cout << " finish\n";
    if (fileIn.peek() != EOF) { std::cout << "Max sequence count reached.\n"; }
    fileIn.close();
  }
  {  // 2. 从索引计算packed偏移
    indices.shrink_to_fit();  // 紧凑化
    std::cout << "sort:\t" << std::flush;
    std::stable_sort(indices.begin(), indices.end(),
      [](const seqIndex& a, const seqIndex& b) {  // 长度递减 稳定排序
        return a.readLength > b.readLength;
      }
    );
    std::cout << "finish\n";
    std::cout << "index:\t" << std::flush;
    size_t baseOffset = 0;  // 基础偏移
    baseOffset += sizeof(uint32_t)
      + indices.size() * sizeof(uint32_t) * 2
      + indices.size() * sizeof(size_t) * 2;
    baseOffset += indices.size() * sizeof(uint32_t) * 128;
    for (uint32_t i = 0; i < indices.size(); i++) {
      indices[i].packedOffsetPack = baseOffset;
      uint32_t byteCount = 1 + (indices[i].readLength + 31) / 32 * 4;
      baseOffset += byteCount * sizeof(uint32_t);
    }
    for (uint32_t i = 0; i < indices.size(); i++) {
      indices[i].packedOffsetSeq = baseOffset;
      baseOffset += indices[i].nameLength + indices[i].readLength + 2;
    }
    std::cout << "finish\n";
  }
  {  // 输出信息
    std::cout << "readed:\t" << indices.size() << "\n";
    std::cout << "whole:\t" << seqCount << "\n";
  }
  return indices;
}

// 按照indices的结果，多线程的写入packed文件
void packData(const std::string& fastaFile, const std::string& packedFile,
const std::vector<seqIndex>& indices) {
  // 精简14转码表 http://bioinfor.imu.edu.cn/raacbook/public/info.html?id=49
  static const std::array<uint32_t, 256> transTable = []() {
    std::array<uint32_t, 256> t = {};
    t['a']= 1; t['c']= 2; t['d']= 3; t['e']= 4; t['f']= 5; t['g']= 6; t['h']= 7;
    t['A']= 1; t['C']= 2; t['D']= 3; t['E']= 4; t['F']= 5; t['G']= 6; t['H']= 7;
    t['i']= 8; t['k']= 9; t['l']=10; t['m']=10; t['n']=11; t['p']=12; t['q']= 4;
    t['I']= 8; t['K']= 9; t['L']=10; t['M']=10; t['N']=11; t['P']=12; t['Q']= 4;
    t['r']=13; t['s']=14; t['t']=15; t['v']= 8; t['w']= 5; t['y']= 5;
    t['R']=13; t['S']=14; t['T']=15; t['V']= 8; t['W']= 5; t['Y']= 5;
    return t;
  }();
  const uint32_t seqCount = indices.size();  // 序列数
  {  // 1. 写入文件头部：序列数、name/read净长度、packed偏移表、fasta偏移表
    std::ofstream fileOut(packedFile, std::ios::binary);
    fileOut.write((char*)&seqCount, sizeof(uint32_t));  // 序列数
    std::vector<uint32_t> nameLengths(seqCount);  // 序列名字净长度
    std::vector<uint32_t> readLengths(seqCount);  // 序列数据净长度
    std::vector<size_t> packOffsets(seqCount);  // 打包偏移
    std::vector<size_t> seqOffsets(seqCount);  // fasta偏移
    for (uint32_t i = 0; i < seqCount; i++) {
      nameLengths[i] = indices[i].nameLength;  // 包含>
      readLengths[i] = indices[i].readLength;
      packOffsets[i] = indices[i].packedOffsetPack;  // 打包偏移
      seqOffsets[i] = indices[i].packedOffsetSeq;  // fasta偏移
    }
    fileOut.write((char*)nameLengths.data(), seqCount * sizeof(uint32_t));
    fileOut.write((char*)readLengths.data(), seqCount * sizeof(uint32_t));
    fileOut.write((char*)packOffsets.data(), seqCount * sizeof(size_t));
    fileOut.write((char*)seqOffsets.data(), seqCount * sizeof(size_t));
  }
  {  // 2. 并发打包和签名，结果直接写入文件
    std::cout << "pack:\t" << std::flush;
    const size_t signBase = sizeof(uint32_t) + seqCount * sizeof(size_t) * 2
      + seqCount * sizeof(uint32_t) * 2;  // 签名的起始偏移
    #pragma omp parallel
    {  // 每个线程有单独的私有变量
      std::ifstream fileIn(fastaFile);  // fasta文件
      std::fstream fileOut(packedFile,
      std::ios::binary | std::ios::in | std::ios::out);  // packed文件
      std::vector<uint32_t> packedBuffer(65536 / 32 * 4 + 1);  // 打包后序列
      std::vector<uint32_t> hashSigns(128);  // 128个哈希签名
      std::string name="", read="", line="";  // 名字、序列、行
      #pragma omp for schedule(dynamic, 1)
      for (uint32_t i = 0; i < seqCount; i++) {  // 遍历序列
        {  // 2.1 从fasta文件读取一条完整记录并解析
          fileIn.seekg(indices[i].fastaOffset);  // 跳转到序列起始位置
          line.resize(indices[i].seqLength);  // 分配内存
          fileIn.read(&line[0], indices[i].seqLength);  // 读数据
          size_t pos = line.find('\n');  // 分割名字和序列
          name = line.substr(0, pos);  // 名字
          read = line.substr(pos + 1);  // 序列
          read.erase(std::remove(read.begin(), read.end(), '\n'), read.end());
        }
        {  // 2.2 打包 + 签名
          const uint32_t length = (uint32_t)read.size();  // 序列长度
          packedBuffer.assign(1 + (length + 31) / 32 * 4, 0);  // 分配内存
          hashSigns.assign(128, 0xFFFFFFFF);  // 重置签名
          uint32_t kmer = 0;  // k-mer 这里是6
          packedBuffer[0] = length;
          for (uint32_t j = 0; j < length; j++) {
            uint32_t pack = transTable[(uint8_t)read[j]];  // 查表
            for (uint32_t e = 0; e < 4; e++) {  // 打包
              packedBuffer[1 + j / 32 * 4 + e] += (pack >> e & 1) << (j % 32);
            }
            kmer <<= 4;
            kmer += pack;
            kmer &= (1 << 25) - 1;  // 截取固定位
            if (j < 6 && j < length - 1) { continue; }  // 前6个不签名
            for (uint32_t k = 0; k < 128; k++) {  // 128个签名
              uint32_t seed = k * 0x9e3779b1 + 0x85ebca6b;  // 黄金比例常数
              uint32_t sign = kmer ^ seed;  // wang's 哈希
              sign = (sign ^ 61) ^ (sign >> 16);
              sign = sign + (sign << 3);
              sign = sign ^ (sign >> 4);
              sign = sign * 0x27d4eb2d;
              sign = sign ^ (sign >> 15);
              sign = (sign ^ (sign >> 8)) * 0x9e3779b1;  // 额外的雪崩效应
              sign = sign ^ (sign >> 14);
              if (sign < hashSigns[k]) { hashSigns[k] = sign; }  // 更新签名
            }
          }
        }
        {  // 2.3 写入packed文件
          fileOut.seekp(signBase + i * 128 * sizeof(uint32_t));  // 写签名
          fileOut.write((char*)hashSigns.data(), 128 * sizeof(uint32_t));
          fileOut.seekp(indices[i].packedOffsetPack);  // 写打包
          uint32_t lengthTemp = packedBuffer.size() * sizeof(uint32_t);
          fileOut.write((char*)packedBuffer.data(), lengthTemp);
          fileOut.seekp(indices[i].packedOffsetSeq);  // 写fasta数据
          line.clear(); line += name; line += "\n"; line += read; line += "\n";
          fileOut.write(line.data(), line.size());
        }
        // 2.4 打印进度
        if (i % (1024 * 1024) == 0) {
          #pragma omp critical
          { std::cout << "." << std::flush; }
        }
      }
    }
    std::cout << " finish\n";  // 进度条结束
  }
}

//--------主函数(仅展示流程)--------//
void makeDB(int argc, char **argv) {
  Timer::Timer timer;
  timer.printStamp();

  // 1. 初始化
  Option option = init(argc, argv);
  // 2. 构建索引
  std::vector<seqIndex> indices = makeIndex(option.fastaFile);
  // 3. 打包数据
  packData(option.fastaFile, option.packedFile, indices);

  timer.printStamp();
  timer.printDuration();
}

}  // namespace MakeDB
#endif  // MAKEDB_HPP
