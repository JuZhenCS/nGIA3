#include <iostream>  // cout
#include "makedb.hpp"  // makeDB
#include "clustering.hpp"  // clustering

int main(int argc, char **argv) {
  if (argc>1 && std::string(argv[1]) == "makedb") {  // 生成数据库
    MakeDB::makeDB(argc-1, argv+1);
  } else if (argc>1 && std::string(argv[1]) == "clustering") {  // 聚类
    Clustering::clustering(argc-1, argv+1);
  } else {  // 打印用法
    std::cout << "Usage:\n" << argv[0] << " makedb/clustering\n";
  }
  return 0;
}

// 今天上午干点啥呢？
// 争取把序列比对部分校对完
// 步骤分解：
// 1. 先读懂test部分，解码序列及之前的内容
// 2. 让test的比对部分，输出比对顺序和结果
// 3. 用python算一下，对应的比对结果，看看是否一致

// 4. 把读数据部分独立出来
// 5. 改写clustering代码，输出比对的顺序和结果
// 6. 校验clustering代码的结果，看看是否一致

// 6. test部分加入阈值判断，然后直接给出比对结果
// 7. clustering加入阈值判断，也直接给出比对结果
// 8. clustering优化左右shift，继续给出比对结果
// 9. test保存结果
// 10. clustering保存结果
