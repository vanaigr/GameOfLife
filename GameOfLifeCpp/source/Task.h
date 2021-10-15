#pragma once
#include <atomic>
#include<thread>
#include<chrono>
#include<cassert>
#include <stdint.h>

#include<iostream>

static size_t counter = 0;

template<class Data>
class Task
{
private:
	
public:
	Data data;
private:
	void(*job)(Data&);
	std::atomic_bool working{false}; //false - waiting, true - working
	std::thread thread;
	size_t index;
	std::chrono::time_point<std::chrono::steady_clock> startTime = std::chrono::steady_clock::now();
public:
	Task(void(*job_)(Data&), Data data_);

	Task(const Task&) = delete;
	Task& operator=(Task const&) = delete;
	~Task() = default;
public:
	void start();
	void waitForResult();
private:
	void task_();
};

template<class Data>
void Task<Data>::task_() {
	while (true) {
		while (working.load() != true) {} // { std::this_thread::sleep_for(std::chrono::microseconds(100)); } //spin
		job(data);
		//std::cout << "work ended " << index << " ms:" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime).count() << std::endl;
		startTime = std::chrono::steady_clock::now();
		working.store(false);
	}
}

template<class Data>
Task<Data>::Task(void(*job_)(Data&), Data data_) : job(job_), data(std::move(data_)), thread{ &Task::task_, this } {
	index = counter++;
}

template<class Data>
void Task<Data>::start() {
	auto expected = false;
	bool isSet = working.compare_exchange_strong(expected, true);
	//std::cout << "start " << index << std::endl;
	startTime = std::chrono::steady_clock::now();
	assert(isSet); //can't do more than 1 task at a time
}

template<class Data>
void Task<Data>::waitForResult() {
	while (working.load() != false) {} //spin
	//std::cout << "spin ended " << index << " ms:" << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime).count() << std::endl;
}