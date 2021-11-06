#pragma once
#include <atomic>
#include<thread>
#include<chrono>
#include<cassert>
#include <stdint.h>

#include<iostream>
#include<mutex>
#include<condition_variable>

static size_t counter = 0;

class Data {

};
template<class Data>
class Task
{
private:
	
public:
	Data data;
private:
	void(*job)(Data&);
	std::mutex startLock;
	std::condition_variable startWork;
	std::atomic_bool workStarted{ false }, workEnded{ false };

	std::thread thread;
	std::chrono::time_point<std::chrono::steady_clock> startTime = std::chrono::steady_clock::now();
public:
	Task(void(*job_)(Data&), Data data_);

	Task(const Task&) = delete;
	Task& operator=(Task const&) = delete;
	~Task() {
		waitForResult();
	}
public:
	void start() noexcept;
	void waitForResult() noexcept;
private:
	void task_() noexcept;
};

template<class Data>
void Task<Data>::task_() noexcept {
	while (true) {
		std::unique_lock<std::mutex> lock{ startLock };
		startWork.wait(lock , [this]() { return workStarted.load() == true; });
		job(data);
		workStarted.store(false);
		workEnded.store(true);
	}
}

template<class Data>
Task<Data>::Task(void(*job_)(Data&), Data data_) : data(std::move(data_)), job(job_), thread{ &Task::task_, this } {}

template<class Data>
void Task<Data>::start() noexcept {
	std::lock_guard<std::mutex> lk{ startLock };
	workEnded.store(false);
	workStarted.store(true);
	startWork.notify_all();
}

template<class Data>
void Task<Data>::waitForResult() noexcept {
	while (workEnded.load() == false) {}

}