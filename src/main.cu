#include <iostream>  // cout
#include "makedb.hpp"  // makeDB

int main(int argc, char **argv) {
  if (argc>1 && std::string(argv[1]) == "makedb") {  // 生成数据库
    MakeDB::makeDB(argc-1, argv+1);
  } else if (argc>1 && std::string(argv[1]) == "clustering") {  // 聚类
    std::cout << "clustering\n";
  } else {  // 打印用法
    std::cout << "Usage:\n" << argv[0] << " makedb/clustering\n";
  }
  return 0;
}

