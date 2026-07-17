#include <exception>
#include <iostream>
#include <string>

#include "makedb.hpp"

int main(int argc, char **argv) {
  try {
    if (argc > 1 && std::string(argv[1]) == "makedb") {
      MakeDB::makeDB(argc - 1, argv + 1);
    } else if (argc > 1 && std::string(argv[1]) == "clustering") {
      std::cout << "clustering\n";
    } else {
      std::cerr << "Usage: " << argv[0] << " makedb -f INPUT -p OUTPUT\n";
      return 2;
    }
  } catch (const std::exception& error) {
    std::cerr << "ngia3: error: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
