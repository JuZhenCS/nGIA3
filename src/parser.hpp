/*
  parser.hpp
  命令行参数解析器 用法如下
  Parser::Parser parser;  // 声明
  parser.add(必要性, 参数名, 参数描述, 默认值);  注册参数
  parser.parse(argc, argv);  解析参数
  getValue<T>(参数名);  读取参数
  没有错误处理，一定不能瞎用！
  2024/05/09 by 鞠震
*/

#ifndef __PARSERHPP__  // 防御声明 pragma once 通用性有问题
#define __PARSERHPP__
#include <iostream>  // cout
#include <vector>  // vector
#include <sstream>  // stringstream
//--------类声明--------//
namespace Parser {  // 命名空间
class Parser {  // 解析器
private:
  struct Parameter {  // 参数
    bool necessary;  // 是否必要
    std::string name;  // 名字
    std::string describe;  // 描述
    std::string value;  // 值
  };
  std::vector<Parameter> parameters;  // 参数集
public:
  Parser() { this->parameters.clear(); }  // 构建
  ~Parser() { this->parameters.clear(); }  // 析构
  void add(const bool necessary, const std::string name,
  const std::string describe, std::string value){  // 添加参数
    this->parameters.push_back({necessary, name, describe, value});
  }
  bool parse(int argc, char **argv) {  // 解析参数
    for (Parameter &parameter : this->parameters) {  // 遍历参数
      for (int i=1; i<argc-1; i++) {  // 找参数
        if (parameter.name == argv[i]) { parameter.value = argv[i+1]; }  // 赋值
      }
      if (parameter.necessary && parameter.value=="") {  // 必要参数未满足
        std::cout << "usage: " << argv[0] << " ... (* is necessary)\n";
        for (const Parameter &parameter : this->parameters) {  // 打印必须
          if (parameter.necessary) {
            std::cout<<"* "<<parameter.name<<"\t"<<parameter.describe<<"\n"; 
          }
        }
        for (const Parameter &parameter : this->parameters) {  // 打印可选
          if (!parameter.necessary) {
            std::cout<<"  "<<parameter.name<<"\t"<<parameter.describe<<"\n";
          }
        }
        return false;
      }
    }
    return true;
  }
  template<typename T>
  T getValue(const std::string &name) { // 读取参数
    std::string content = "";  // 内容
    for (const Parameter &parameter : this->parameters) {  // 取内容
      if (parameter.name == name) { content = parameter.value; }
    }
    T value;  // 转换后的值
    std::istringstream iss(content);  // 输入流
    iss >> value;  // 转换
    return value;  // 返回值
  }
};
} // namespace Parser
#endif // __PARSERHPP__
