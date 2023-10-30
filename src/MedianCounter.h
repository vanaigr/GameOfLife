#pragma once
#include<deque>
#include<iostream>

class MedianCounter {
private:
	double accumulator = 0;
	unsigned int count = 0;
public:
	MedianCounter() = default;

	void add(double val) {
		accumulator += val;
		count++;
	}

	double median() const {
		return accumulator / count;
	}
};

class UMedianCounter {
private:
	const size_t size;
	std::deque<double> accumulator;
public:
	UMedianCounter(size_t size_) : size{ size_ }, accumulator{} {};

	void add(double val) {
		if (accumulator.size() == size) accumulator.pop_front();
		accumulator.push_back(val);
	}

	double median() const {
		double value = 0;
		for (size_t i = 0; i < accumulator.size(); i++) {
			value += accumulator[i];
		}
		return value / accumulator.size();
	}

	double min() const {
		if (accumulator.size() == 0) {
			return 0.0;
		}
		else {
			double min = accumulator[0];
			for (size_t i = 1; i < accumulator.size(); i++) {
				if(min > accumulator[i]) min = accumulator[i];
			}
			return min;
		}
	}

	double max() const {
		if (accumulator.size() == 0) {
			return 0.0;
		}
		else {
			double max = accumulator[0];
			for (size_t i = 1; i < accumulator.size(); i++) {
				if (max < accumulator[i]) max = accumulator[i];
			}
			return max;
		}
	}

	void printall(std::ostream &stream) const {
		for (auto el : accumulator) {
			stream << el << std::endl;
		}
	}
};


