#pragma once

#include<iostream>
#include <stdint.h>
#include<memory>
#include "Vector.h"
#include <thread>
#include <atomic>
#include "Task.h"
#include <atomic>
#include <vector>
#include"MedianCounter.h"
#include"Misc.h"
#include<unordered_set>
#include<functional>

using FieldCell = bool;

namespace fieldCell {
    static constexpr FieldCell cellAlive = true;
    static constexpr FieldCell cellDead = false;

    inline FieldCell nextGeneration(FieldCell const cell, uint8_t const aliveNeighboursCount) {
        if(!cell && aliveNeighboursCount == 3) return true;
        else if(cell && aliveNeighboursCount != 2 && aliveNeighboursCount != 3) return false;
        else return cell;
    }
    inline constexpr char const *asString(const FieldCell cell) {
        return cell ? "alive" : "dead";
    }
}

struct FieldModification {
    int64_t startIndex_int, size_int;
    const uint32_t *data;
};

struct FieldOutput { //TODO: remove
    virtual void write(FieldModification fm) {
        this->batched()->write(fm);
    }

    virtual std::unique_ptr<FieldOutput> batched() const = 0;
    virtual std::unique_ptr<FieldOutput> tryBatched() const = 0;
    virtual ~FieldOutput() = default;
}; //must be used as output only in one thread


struct Cell {
    FieldCell cell;
    int64_t index;
};

class Field final {
public:
    struct FieldPimpl;
    struct GridData;
private:
    std::unique_ptr<FieldPimpl> gridPimpl;
    bool isStopped;
    std::unique_ptr<FieldOutput> const current_output;
    std::unique_ptr<FieldOutput> const buffer_output;

    int numberOfTasks;
    std::unique_ptr<std::unique_ptr<Task<GridData>>[/*numberOfTasks*/]> gridTasks;
    std::atomic_bool interrupt_flag;
    std::unordered_set<int64_t> indecesToBrokenCells;
public:
    Field(
        const int64_t gridWidth, const int64_t gridHeight, const int numberOfTasks_,
        std::function<std::unique_ptr<FieldOutput>()> current_outputs,
        std::function<std::unique_ptr<FieldOutput>()> buffer_outputs
    );
    ~Field();

    Field(Field const&) = delete;
    Field& operator=(Field const&) = delete;
public:
    bool tryFinishGeneration();
    void startCurGeneration();
    void startNewGeneration();

    void fill(const FieldCell cell);

    FieldCell cellAtIndex(const int64_t index) const;

    void setCellAtIndex(const int64_t index, FieldCell cell);
    void setCellAtCoord(const vec2l& coord, FieldCell cell);
    void setCells(Cell const *const cells, size_t const count);

    FieldCell cellAtCoord(const vec2l& coord) const;
    FieldCell cellAtCoord(const int64_t column, const int64_t row) const;

    vec2l indexAsCoord(const int64_t index) const;

    int64_t coordAsIndex(const vec2l& coord) const;
    int64_t coordAsIndex(const int64_t column, const int64_t row) const;

    int64_t normalizeIndex(const int64_t index) const;
    vec2l normalizeCoord(const vec2l& coord) const;

    int64_t width() const;
    int64_t height() const;
    int64_t size() const;

    int64_t size_bytes() const;
    int64_t width_actual() const;
    uint32_t *rawData() const;

    void halt();
private:
    void waitForGridTasks();
    void deployGridTasks();
};

inline void Field::setCellAtCoord(const vec2l& coord, FieldCell cell) {
    setCellAtIndex(coordAsIndex(coord), cell);
}

inline FieldCell Field::cellAtCoord(const vec2l& coord) const {
    return cellAtIndex(coordAsIndex(coord));
}

inline FieldCell Field::cellAtCoord(const int64_t column, const int64_t row) const {
    return cellAtCoord(vec2l(column, row));
}

inline vec2l Field::indexAsCoord(const int64_t index) const {
    const auto index_n = normalizeIndex(index);
    const auto x = index_n % width();
    const auto y = index_n / width();
    return vec2l(x, y);
}

inline int64_t Field::coordAsIndex(const vec2l& coord) const {
    const auto coord_n = normalizeCoord(coord);
    return (int64_t)coord_n.x + (int64_t)coord_n.y * (int64_t)width();
}

inline int64_t Field::coordAsIndex(const int64_t column, const int64_t row) const {
    return coordAsIndex(vec2l(column, row));
}

inline int64_t Field::normalizeIndex(const int64_t index) const {
    return misc::mod(index, size());
}
inline vec2l Field::normalizeCoord(const vec2l& coord) const {
    return vec2l(misc::mod(coord.x, width()), misc::mod(coord.y, height()));
}

