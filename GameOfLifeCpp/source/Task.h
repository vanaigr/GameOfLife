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
	volatile bool continueThread { true };

	std::thread thread;
	std::chrono::time_point<std::chrono::steady_clock> startTime = std::chrono::steady_clock::now();
public:
	template<class... DataArgs>
	Task(void(*job_)(Data&), DataArgs&&... args) : data(std::forward<DataArgs>(args)...), job(job_), thread{ &Task::task_, this } {}

	Task(const Task&) = delete;
	Task& operator=(Task const&) = delete;
	~Task() noexcept {
		{
			std::lock_guard<std::mutex> lk{ startLock };
			continueThread = false;
			startWork.notify_one();
		}
		thread.join();
	}
public:
	void start() noexcept;
	void waitForResult() noexcept;
private:
	void task_() noexcept;
};

template<class Data>
void Task<Data>::task_() noexcept {
	bool exit_ = continueThread;
	while (exit_) {
		std::unique_lock<std::mutex> lock{ startLock };
		startWork.wait(lock , [this, &exit_]() { return workStarted.load() == true || (exit_ = !continueThread); });
		if (exit_) return;
		job(data);
		workStarted.store(false);
		workEnded.store(true);
	}
}

template<class Data>
void Task<Data>::start() noexcept {
	std::lock_guard<std::mutex> lk{ startLock };
	workEnded.store(false);
	workStarted.store(true);
	startWork.notify_one();
}

template<class Data>
void Task<Data>::waitForResult() noexcept {
	while (workEnded.load() == false) {}

}