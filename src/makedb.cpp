/*
makedb.cpp
输出文件数据内容:
  uint32_t         序列数         1
  uint32_t[seqCount] nameLengths   序列名长度
  uint32_t[seqCount] readLengths   序列长度
  vector<uint64_t> packed数据偏移  seqCount
  vector<uint64_t> fasta数据偏移   seqCount
  MinHash签名      uint32_t[128] * seqCount
  packed数据 (readsCount条记录)
    uint32_t 数据长度 1
    uint32_t 压缩数据 (数据长度+31)/32*4
  fasta数据 (readsCount条记录)
    string 序列名 1
    string 序列 1
长度<0xFFFF 数量<=0x7FFFFFFF
用法:
ngia3 makedb -f fasta文件 -p packed文件
2026-07-09 by 鞠震
*/

// 序列长 < 0xFFFF     比对时line  < 2048
// 序列数 < 0x7FFFFFFF 比对时线程数 < int范围


#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

#include "makedb.hpp"
#include "timer.hpp"

namespace MakeDB {  // 命名空间
//--------------------------------  数据结构  --------------------------------//

namespace {
struct Option {  // 输入选项
  std::string fastaFile;  // fasta文件
  std::string packedFile; // packed文件
};

struct SequenceIndex {  // 序列索引
  uint64_t fastaOffset;  // fasta文件偏移
  uint64_t packedOffsetPack;  // packed文件中 打包偏移
  uint64_t packedOffsetSeq;  // packed文件中 序列偏移
  uint32_t nameLength;  // 序列名字长度
  uint32_t readLength;  // 序列数据长度
  uint64_t seqLength;  // 序列记录长度 名字+数据+换行
};

inline size_t lineContentLength(const std::string& line) {
  const size_t carriageReturn =
    (!line.empty() && line.back() == '\r') ? 1 : 0;
  return line.size() - carriageReturn;
}

constexpr uint32_t MAX_SEQUENCE_LENGTH = 0xFFFF;
constexpr uint32_t MAX_SEQUENCE_COUNT = 0x7FFFFFFF;
constexpr size_t SIGNATURE_COUNT = 128;

inline const std::array<uint32_t, SIGNATURE_COUNT>& hashSeeds() {
  static const std::array<uint32_t, SIGNATURE_COUNT> seeds = []() {
    std::array<uint32_t, SIGNATURE_COUNT> result{};
    for (size_t index = 0; index < result.size(); ++index) {
      result[index] =
        static_cast<uint32_t>(index) * 0x9e3779b1U + 0x85ebca6bU;
    }
    return result;
  }();
  return seeds;
}

inline uint32_t wangHash(uint32_t value, uint32_t seed) {
  value ^= seed;
  value = (value ^ 61U) ^ (value >> 16);
  value += value << 3;
  value ^= value >> 4;
  value *= 0x27d4eb2dU;
  value ^= value >> 15;
  value = (value ^ (value >> 8)) * 0x9e3779b1U;
  return value ^ (value >> 14);
}


//--------------------------------  功能函数  --------------------------------//
// init初始化
inline Option init(int argc, char** argv) {
  Option option{};
  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    if (argument != "-f" && argument != "-p") {
      throw std::invalid_argument("unknown makedb argument: " + argument);
    }
    if (index + 1 >= argc) {
      throw std::invalid_argument("missing value for " + argument);
    }
    const std::string value = argv[++index];
    if (argument == "-f") {
      if (!option.fastaFile.empty()) {
        throw std::invalid_argument("duplicate -f argument");
      }
      option.fastaFile = value;
    } else {
      if (!option.packedFile.empty()) {
        throw std::invalid_argument("duplicate -p argument");
      }
      option.packedFile = value;
    }
  }
  if (option.fastaFile.empty() || option.packedFile.empty()) {
    throw std::invalid_argument("makedb requires -f INPUT and -p OUTPUT");
  }
  if (std::filesystem::exists(option.packedFile) &&
      std::filesystem::equivalent(option.fastaFile, option.packedFile)) {
    throw std::invalid_argument(
      "input FASTA and packed output must be different files");
  }
  bool checkPassed = true;  // 校验是否通过
  {  // 校验参数
    // 1. 文件存在
    std::ifstream fastaFile(option.fastaFile, std::ios::binary);  // fasta文件
    if (!fastaFile.good()) {  // 文件不存在
      std::cout << option.fastaFile << " does not exist.\n";
      checkPassed = false;
    } else {
      char previousChar = '\0';
      char currentChar = '\0';
      char lastChar = '\0';
      bool bareCarriageReturn = false;
      while (fastaFile.get(currentChar)) {
        if (previousChar == '\r' && currentChar != '\n') {
          bareCarriageReturn = true;
        }
        previousChar = currentChar;
        lastChar = currentChar;
      }
      if (previousChar == '\r') {
        bareCarriageReturn = true;
      }
      if (bareCarriageReturn) {
        std::cout << "File contains a bare \\r character.\n";
        checkPassed = false;
      }
      if (lastChar != '\0' && lastChar != '\n') {
        std::cout << "File must end with \\n.\n";
        checkPassed = false;
      }
    }
  }
  if (!checkPassed) { throw std::runtime_error("input validation failed"); }
  std::cout << "fasta:\t" << option.fastaFile << "\n";  // 输入数据
  std::cout << "packed:\t" << option.packedFile << "\n";  // 打包输出
  std::cout << "thread:\t" << omp_get_max_threads() << "\n";  // 可用线程数
  return option;
}

// makeIndex 第一遍读文件，构建索引
inline std::vector<SequenceIndex> makeIndex(const std::string& fastaFile) {
  std::vector<SequenceIndex> indices;  // 序列索引
  uint32_t seqCount = 0;  // 扫描的序列数（含超长丢弃的）
  {  // 1. 从文件读序列索引
    std::cout << "read:\t." << std::flush;  // 进度条
    std::ifstream fileIn(fastaFile, std::ios::binary);  // fasta文件
    std::string line;  // 一行数据
    uint64_t fastaOffset = 0;  // fasta偏移
    uint64_t nameLength = 0;  // 序列名字长度
    uint64_t readLength = 0;  // 序列数据长度
    uint64_t seqLength = 0;  // 序列长度
    while (fileIn.peek() != EOF && seqCount < MAX_SEQUENCE_COUNT) {
      fastaOffset = static_cast<uint64_t>(fileIn.tellg());  // fasta偏移
      getline(fileIn, line);  // 读 >name 行
      if (line.empty() || line.front() != '>') {
        throw std::runtime_error(
          "invalid FASTA header at record " + std::to_string(seqCount + 1));
      }
      nameLength = lineContentLength(line);  // 序列名字长度
      seqLength = line.size() + 1;  // 包括换行符
      readLength = 0;  // 序列数据长度
      while (fileIn.peek() != EOF && fileIn.peek() != '>') {
        getline(fileIn, line);  // 序列数据
        readLength += lineContentLength(line);  // 序列数据长度
        seqLength += line.size() + 1;  // 包括换行符
      }
      seqCount += 1;  // 文件中的序列总数
      if (readLength < MAX_SEQUENCE_LENGTH) {
        if (nameLength > std::numeric_limits<uint32_t>::max()) {
          throw std::runtime_error("FASTA header exceeds uint32_t length");
        }
        indices.push_back({
          fastaOffset, 0, 0, static_cast<uint32_t>(nameLength),
          static_cast<uint32_t>(readLength), seqLength
        });
      }
      if (seqCount % (1024 * 1024) == 0) { std::cout << "." << std::flush; }
    }
    std::cout << " finish\n";
    if (fileIn.peek() != EOF) { std::cout << "Max sequence count reached.\n"; }
  }
  {  // 2. 从索引计算packed偏移
    std::cout << "sort:\t" << std::flush;
    std::stable_sort(indices.begin(), indices.end(),
      [](const SequenceIndex& a, const SequenceIndex& b) {
        return a.readLength > b.readLength;
      }
    );
    std::cout << "finish\n";
    std::cout << "index:\t" << std::flush;
    uint64_t baseOffset = sizeof(uint32_t);
    baseOffset += indices.size() * (sizeof(uint32_t) * 2 + sizeof(uint64_t) * 2);
    baseOffset += indices.size() * sizeof(uint32_t) * SIGNATURE_COUNT;
    for (size_t i = 0; i < indices.size(); ++i) {
      indices[i].packedOffsetPack = baseOffset;
      uint32_t byteCount = 1 + (indices[i].readLength + 31) / 32 * 4;
      baseOffset += byteCount * sizeof(uint32_t);
    }
    for (size_t i = 0; i < indices.size(); ++i) {
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
inline void packData(const std::string& fastaFile, const std::string& packedFile,
                     const std::vector<SequenceIndex>& indices) {
  const uint16_t byteOrderMarker = 1;
  if (*reinterpret_cast<const uint8_t*>(&byteOrderMarker) != 1) {
    throw std::runtime_error(
      "the nGIA3 packed format currently requires a little-endian host");
  }
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
  if (indices.size() > std::numeric_limits<uint32_t>::max()) {
    throw std::runtime_error("too many retained sequences");
  }
  const uint32_t seqCount = static_cast<uint32_t>(indices.size());
  std::fstream fileOut(
    packedFile,
    std::ios::binary | std::ios::in | std::ios::out | std::ios::trunc);
  if (!fileOut) { throw std::runtime_error("cannot open packed output"); }
  fileOut.write(reinterpret_cast<const char*>(&seqCount), sizeof(seqCount));
  std::vector<uint32_t> nameLengths(seqCount);
  std::vector<uint32_t> readLengths(seqCount);
  std::vector<uint64_t> packOffsets(seqCount);
  std::vector<uint64_t> seqOffsets(seqCount);
  for (uint32_t i = 0; i < seqCount; ++i) {
    nameLengths[i] = indices[i].nameLength;
    readLengths[i] = indices[i].readLength;
    packOffsets[i] = indices[i].packedOffsetPack;
    seqOffsets[i] = indices[i].packedOffsetSeq;
  }
  fileOut.write(reinterpret_cast<const char*>(nameLengths.data()),
                seqCount * sizeof(uint32_t));
  fileOut.write(reinterpret_cast<const char*>(readLengths.data()),
                seqCount * sizeof(uint32_t));
  fileOut.write(reinterpret_cast<const char*>(packOffsets.data()),
                seqCount * sizeof(uint64_t));
  fileOut.write(reinterpret_cast<const char*>(seqOffsets.data()),
                seqCount * sizeof(uint64_t));
  if (!fileOut) { throw std::runtime_error("failed to write packed header"); }

  std::cout << "pack:\t" << std::flush;
  const uint64_t groupBase = sizeof(uint32_t) +
    seqCount * (sizeof(uint32_t) * 2 + sizeof(uint64_t) * 2);
  std::atomic_bool ioFailed{false};
  #pragma omp parallel shared(fileOut, ioFailed)

    {  // 每个线程有单独的私有变量
      std::ifstream fileIn(fastaFile, std::ios::binary);  // fasta文件
      if (!fileIn) { ioFailed.store(true, std::memory_order_relaxed); }
      std::vector<uint32_t> packedBuffer;
      packedBuffer.reserve(1 + (MAX_SEQUENCE_LENGTH + 31) / 32 * 4);
      std::array<uint32_t, SIGNATURE_COUNT> hashSigns{};
      const auto& seeds = hashSeeds();
      std::string name, read, line;
      #pragma omp for schedule(dynamic, 1)
      for (uint32_t i = 0; i < seqCount; ++i) {
        if (ioFailed.load(std::memory_order_relaxed)) {
          continue;
        }
        {  // 2.1 从fasta文件读取一条完整记录并解析
          fileIn.clear();
          fileIn.seekg(static_cast<std::streamoff>(indices[i].fastaOffset));
          const std::streamsize recordSize =
            static_cast<std::streamsize>(indices[i].seqLength);
          line.assign(static_cast<size_t>(recordSize), '\0');
          fileIn.read(line.data(), recordSize);
          if (fileIn.gcount() != recordSize) {
            ioFailed.store(true, std::memory_order_relaxed);
            continue;
          }
          const size_t pos = line.find('\n');
          name = line.substr(0, pos);
          if (!name.empty() && name.back() == '\r') { name.pop_back(); }
          read = line.substr(pos + 1);
          read.erase(std::remove(read.begin(), read.end(), '\n'), read.end());
          read.erase(std::remove(read.begin(), read.end(), '\r'), read.end());
        }
        {  // 2.2 打包 + 签名
          const uint32_t length = static_cast<uint32_t>(read.size());
          packedBuffer.assign(1 + (length + 31) / 32 * 4, 0);  // 分配内存
          hashSigns.fill(std::numeric_limits<uint32_t>::max());
          uint32_t kmer = 0;  // k-mer 这里是6
          packedBuffer[0] = length;
          for (uint32_t j = 0; j < length; ++j) {
            const uint32_t pack = transTable[static_cast<uint8_t>(read[j])];
            for (uint32_t e = 0; e < 4; e++) {  // 打包
              packedBuffer[1 + j / 32 * 4 + e] += (pack >> e & 1) << (j % 32);
            }
            kmer <<= 4;
            kmer += pack;
            kmer &= (1U << 25) - 1U;
            if (j < 6 && j < length - 1) { continue; }  // 前6个不签名
            for (size_t k = 0; k < SIGNATURE_COUNT; ++k) {
              const uint32_t sign = wangHash(kmer, seeds[k]);
              if (sign < hashSigns[k]) { hashSigns[k] = sign; }
            }
          }
        }
          line.clear();
        {  // 2.3 写入packed文件
          line.reserve(name.size() + read.size() + 2);
          line.append(name).push_back('\n');
          line.append(read).push_back('\n');
          #pragma omp critical(packed_file_write)
          {
            fileOut.seekp(static_cast<std::streamoff>(
              groupBase + i * SIGNATURE_COUNT * sizeof(uint32_t)));
            fileOut.write(reinterpret_cast<const char*>(hashSigns.data()),
                          SIGNATURE_COUNT * sizeof(uint32_t));
            fileOut.seekp(static_cast<std::streamoff>(
              indices[i].packedOffsetPack));
            const std::streamsize packedBytes = static_cast<std::streamsize>(
              packedBuffer.size() * sizeof(uint32_t));
            fileOut.write(reinterpret_cast<const char*>(packedBuffer.data()),
                          packedBytes);
            fileOut.seekp(static_cast<std::streamoff>(
              indices[i].packedOffsetSeq));
            fileOut.write(line.data(), static_cast<std::streamsize>(line.size()));
            if (!fileOut) { ioFailed.store(true, std::memory_order_relaxed); }
          }
        }
        // 2.4 打印进度
        if (i % (1024 * 1024) == 0) {
          #pragma omp critical(progress_output)
          { std::cout << "." << std::flush; }
        }
      }
    }
  fileOut.flush();
  if (ioFailed.load(std::memory_order_relaxed) || !fileOut) {
    throw std::runtime_error("failed while reading FASTA or writing packed data");
  }
  std::cout << " finish\n";
}


}  // namespace
//--------主函数(仅展示流程)--------//
void makeDB(int argc, char** argv) {
  Timer::Timer timer;
  timer.printStamp();

  // 1. 初始化
  Option option = init(argc, argv);
  // 2. 构建索引
  std::vector<SequenceIndex> indices = makeIndex(option.fastaFile);
  // 3. 打包数据
  packData(option.fastaFile, option.packedFile, indices);

  timer.printStamp();
  timer.printDuration();
}

}  // namespace MakeDB
