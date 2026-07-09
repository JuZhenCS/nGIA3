/*
  timer.hpp
  chrono库实现的计时器 支持时间戳 中断计时 用法如下
  Timer::Timer timer;  声明
  timer.printStamp();  打印时间戳
  timer.start();  开始计时
  timer.pause();  暂停计时
  timer.resume();  恢复计时
  timer.printDuration();  打印耗时
  2025/05/08 by 鞠震
*/

#ifndef __TIMERHPP__  // 防御声明 pragma once 通用性有问题
#define __TIMERHPP__
#include <iostream>  // cout
#include <chrono>  // system_clock
//--------类声明--------//
namespace Timer {  // 命名空间
class Timer {  // 计时器
private:
  std::chrono::steady_clock::time_point startPoint;  // 开始时间
  std::chrono::nanoseconds duration;  // 耗时 纳秒
  bool running;  // 正在计时
public:
  Timer() {  // 构造
    this->startPoint = std::chrono::steady_clock::now();
    this->duration =  std::chrono::nanoseconds(0);
    this->running = true;
  }
  ~Timer() {}  // 析构
  void printStamp() {  // 打印时间戳
    auto stamp = std::chrono::system_clock::now();
    auto formatedStamp = std::chrono::system_clock::to_time_t(stamp);
    std::cout << ctime(&formatedStamp);  // 自带换行
  }
  void start() {  // 开始计时
    this->startPoint = std::chrono::steady_clock::now();
    this->duration = std::chrono::nanoseconds(0);
    this->running = true;
  }
  void pause() {  // 暂停计时
    auto endPoint = std::chrono::steady_clock::now();
    this->duration += endPoint - this->startPoint;
    this->running = false;
  }
  void resume() {  // 恢复计时
    this->startPoint = std::chrono::steady_clock::now();
    this->running = true;
  }
  void printDuration() {  // 打印耗时
    auto endPoint = std::chrono::steady_clock::now();
    auto elapse = this->duration;
    if (running) {elapse += endPoint - this->startPoint;}
    auto hh = std::chrono::duration_cast<std::chrono::hours>(elapse).count();
    auto mm = std::chrono::duration_cast<std::chrono::minutes>(elapse).count();
    auto ss = std::chrono::duration_cast<std::chrono::seconds>(elapse).count();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(elapse).count();
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(elapse).count();
    auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(elapse).count();
    std::cout << "Consuming:\t";  // 计时没有四舍五入
    std::cout << hh << " h\t";
    std::cout << mm%60 << " m\t";
    std::cout << ss%60 << " s\t";
    std::cout << ms%1000 << " ms\t";
    std::cout << us%1000 << " us\t";
    std::cout << ns%1000 << " ns\n";
  }
};
}  // namespace Timer
#endif  // __TIMERHPP__
