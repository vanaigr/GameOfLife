#include"glew.h"
#include<GLFW/glfw3.h>

#include"Misc.h"
#include"Grid.h"
#include<vector>

#include<cassert> 

#include"Timer.h"
#include"AutoTimer.h"
#include"MedianCounter.h"

#include<nmmintrin.h> 

#include<algorithm>

using Cells = uint32_t;
static constexpr auto cellsBatchSize = 4;
static constexpr auto cellsBatchLength = 32;

struct Field::FieldPimpl {
    using index_t = int32_t;

    using BufferType = bool;
    static constexpr BufferType bufCur = false;
    static constexpr BufferType bufNext = true;

    Cells* buffer;

    int32_t width;
    int32_t height;
    int32_t rowLength;
    bool edgeCellsOptimization;
    bool buffersSwapped;


public:
    FieldPimpl(const int32_t gridWidth, const int32_t gridHeight) {
        width = gridWidth;
        height = gridHeight;
        rowLength = misc::intDivCeil(width, cellsBatchLength);
        edgeCellsOptimization = rowLength * cellsBatchLength - width >= 2; 

        auto const bufferLen = bufferLength();
        auto const paddingLen = bufferPaddingLength();
        buffer = new Cells[bufferLen*2]{};
        buffersSwapped = false;
    }
    ~FieldPimpl() = default;

    void fixField(BufferType const type = bufCur) {
        auto const paddingLen = bufferPaddingLength();
        auto const gridLen = gridLength();
        auto const buffer = getBuffer(type) + paddingLen;
        auto const rowSize = rowLength * cellsBatchSize;

        auto const firstEmptyRowCellBatch = width / cellsBatchLength;
        auto const firstEmptyRowCellInBatch = width % cellsBatchLength;

        //order is important as if there can be only one row,
        //in which case this algorithm *should* still work:
        //side neightbours algorithm requires one extra row ahead, which 
        //we get as the end padding IS the next row. then we fix 
        //the side neightbours for start and end padding and copy 
        //start padding (not from -1 but from 0, as -1 can be outside the grid)

        auto const startPaddingRow = buffer - rowLength;
        auto const endPaddingRow = buffer + gridLen;

        //copy end padding row
        std::memcpy(endPaddingRow, buffer, rowSize);

        //copy left/right neighbours to other side
        if(edgeCellsOptimization) {
            for(int32_t row = 0; row != height; row++) {
                auto const rowOffset = row * rowLength;

                auto &cells = *(buffer + rowOffset + firstEmptyRowCellBatch);

                auto const firstCell = *(buffer + rowOffset) & 1;
                auto const lastCell = (*(buffer + rowOffset + rowLength + rowLength-1) >> (firstEmptyRowCellInBatch-1)) & 1;

                cells = (cells & ~(~0 << (firstEmptyRowCellInBatch)))
                    | (firstCell << firstEmptyRowCellInBatch)
                    | (lastCell << (cellsBatchLength-1));
            }

            //copy recalculated neighbours to padding.
            *(startPaddingRow-1) = *(buffer + gridLen - rowLength - 1);
            *(endPaddingRow + rowLength - 1) = *(buffer + rowLength - 1);
        }

        //copy start padding row
        std::memcpy(startPaddingRow, buffer + gridLen - rowLength, rowSize);
    }

    void swapBuffers() { buffersSwapped = !buffersSwapped; }

    void fill(FieldCell const cell, BufferType const type = bufCur) {
        auto grid = getBuffer(type);
        const auto val = ~0u * cell;
        std::fill(&grid[0], &grid[0] + bufferLength(), val);
    }

    uint8_t cellAt_grid(int32_t const index, BufferType const type = bufCur) const {
        auto grid = getBuffer(type);
        const auto row = misc::intDivFloor(index, width);
        const auto col = misc::mod(index, width);

        const auto col_int = col / cellsBatchLength;
        const auto shift = col % cellsBatchLength;

        return (grid[bufferPaddingLength() + row * rowLength + col_int] >> shift) & 0b1;
    }

    void setCellAt(int32_t const index, FieldCell const cell, BufferType const type = bufCur) {
        auto grid = getBuffer(type);

        const auto row = misc::intDivFloor(index, int32_t(width));
        const auto col = misc::mod(index, width);

        const auto col_int = col / cellsBatchLength;
        const auto shift = col % cellsBatchLength;

        auto& cur{ grid[bufferPaddingLength() + row * rowLength + col_int] };

        cur = (cur & ~(0b1u << shift)) | (static_cast<uint32_t>(cell) << shift);
    }

    uint32_t& getCellsActual_int(int32_t const index_actual_int, BufferType const type = bufCur) const {
        return getBuffer(type)[index_actual_int + bufferPaddingLength()];
    }

    uint32_t cellI2BatchI(const uint32_t index) const {
        const auto row = misc::intDivFloor(index, width);
        const auto col = misc::mod(index, width);

        const auto col_int = col / cellsBatchLength;

        return row * rowLength + col_int;
    }

    Cells *getBuffer(BufferType const type) const {
        auto const offset = (type == bufNext) ^ buffersSwapped ? bufferLength() : 0;
        return buffer + offset;
    }
    uint32_t gridLength() const {
        return height * rowLength;
    }
    uint32_t bufferPaddingLength() const {
        return rowLength + 1/*
            extra row before/after the grid repeating the opposite row 
            and 1 cell on each side for the first/last cell's neighbour,
            padding after can be 1 shorter but not when width % batchSize, 
            so it is +1 always for consistency
        */;
    }
    uint32_t bufferLength() const {
        return gridLength() + 2*bufferPaddingLength();
    }
};

struct Field::GridData {
private: static const uint32_t samples = 100;
public:
    UMedianCounter gridUpdate{ samples }, bufferSend{ samples };
    uint32_t task__iteration;
    uint32_t task__index;

    std::unique_ptr<FieldPimpl>& grid;
    std::atomic_bool& interrupt_flag;
    uint32_t startBatch;
    uint32_t endBatch;
    std::unique_ptr<FieldOutput> const buffer_output;

    void generationUpdated() {
        task__iteration++;
        if (task__iteration % (samples + task__index) == 0)
            std::cout << "grid task " << task__index << ": "
            << gridUpdate.median() << ' '
            << bufferSend.median() << std::endl;
    }

public:
    GridData(
        uint32_t index_,
        std::unique_ptr<FieldPimpl>& grid_,
        std::atomic_bool& interrupt_flag_,
        uint32_t startBatch_,
        uint32_t endBatch_,
        std::unique_ptr<FieldOutput> &&output_
    ) :
        task__iteration{ 0 },
        task__index(index_),
        grid(grid_),
        interrupt_flag(interrupt_flag_),
        startBatch(startBatch_),
        endBatch  (endBatch_),
        buffer_output{ std::move(output_) }
    {}
};

struct Remainder {
    uint16_t cellsCols;
    uint8_t curCell;
};
//computes new generation for 32 cells:
//one from previous remainder and 31 cells at *base.
//also computes remainder for next iteration (for the cast cell of *base)
static uint32_t newGenerationBatched(
    Remainder const previousRemainder,
    Cells *base,
    int32_t const rowLength,
    Remainder &currentRemainder_out
) {
    // [0, 16] [1, 17] ... [15, 31], where for cells x, y: [x, y] means 0b000y'000x
    auto const unpackCellsAs4Bits = [](const uint32_t number) -> __m128i {
        auto const cellPosForByteMask = _mm_setr_epi8(
            1, 2, 4, 8, 16, 32, 64, 128u,
            1, 2, 4, 8, 16, 32, 64, 128u
        );

        auto const numberReg = _mm_cvtsi32_si128(number);
        auto const numberDupBytes = _mm_unpacklo_epi8(numberReg, numberReg);

        auto const numQuadrupleLowBytes = _mm_shufflelo_epi16(numberDupBytes, 0b01'01'00'00);
        auto const numberLowHalf = _mm_shuffle_epi32(numQuadrupleLowBytes, 0b01'01'00'00); 
        auto const isLowCell = _mm_cmpeq_epi8(_mm_and_si128(numberLowHalf, cellPosForByteMask), cellPosForByteMask);

        auto const numQuadrupleHighBytes = _mm_shufflelo_epi16(numberDupBytes, 0b11'11'10'10);
        auto const numberHighHalf = _mm_shuffle_epi32(numQuadrupleHighBytes, 0b01'01'00'00); 
        auto const isHighCell = _mm_cmpeq_epi8(_mm_and_si128(numberHighHalf, cellPosForByteMask), cellPosForByteMask);

        return _mm_sub_epi8(_mm_and_si128(isHighCell, _mm_set1_epi8(0b0001'0000)), isLowCell)/*
            low is either:
                true (-1), then  high&16 - -1  =>  high&16 | 1
                false (0), then  high&16 -  0  =>  high&16
            which is eqivalent of  high&16 | low&1
        */;
    };

    auto const topBatch = unpackCellsAs4Bits(*(base - rowLength));
    auto const curBatch = unpackCellsAs4Bits(*base);
    auto const botBatch = unpackCellsAs4Bits(*(base + rowLength));

    auto const verticalSum = _mm_add_epi8(topBatch, _mm_add_epi8(curBatch, botBatch));

    currentRemainder_out.curCell = (_mm_extract_epi16(curBatch, 7) >> 12); //_epi8 is sse 4.1
    currentRemainder_out.cellsCols = (_mm_extract_epi16(verticalSum, 7) >> 4) & 0b1111'00001111;

    auto const first16Mask = _mm_set1_epi8(0b00001111u); //unnecessary if _slli_epi8 existed
    auto const curRowCentered_carry = _mm_slli_epi16(_mm_and_si128(curBatch, first16Mask), 4);
    auto const verticalSum_carry = _mm_slli_epi16(_mm_and_si128(verticalSum, first16Mask), 4);

    //using _or instead of _add because .curCell is in lower 4 bits and carry is in higher 4
    auto const curRowCentered = _mm_or_si128(
        _mm_alignr_epi8(curBatch, curRowCentered_carry, 15), 
        _mm_cvtsi32_si128(previousRemainder.curCell)
    );
    auto const cells3by3 = _mm_add_epi8(
        _mm_add_epi8(
            _mm_add_epi8(
                _mm_alignr_epi8(verticalSum, verticalSum_carry, 14),
                _mm_alignr_epi8(verticalSum, verticalSum_carry, 15)
            ),
            verticalSum // ~ alignr(..., 16)
        ),
        _mm_cvtsi32_si128(previousRemainder.cellsCols + (previousRemainder.cellsCols >> 8))
    );

    auto const cellsNeighboursAlive = _mm_sub_epi8(cells3by3, curRowCentered);
    auto const cells = _mm_or_si128(cellsNeighboursAlive, curRowCentered);
    static_assert(
        (2 | 1) == 3 && (3 | 1) == 3 && (3 | 0) == 3 && true,
        R"(must be true:
            1, 2) (2 or 3 cells) | (alive cell) == 3
            3) (3 cells) | (dead  cell) == 3
            4) other combitations != 3
        )"
    );

    auto const mask_lower = _mm_set1_epi8(0b1111u);
    auto const three = _mm_set1_epi8(3u);

    uint32_t const lower16 = _mm_movemask_epi8(_mm_cmpeq_epi8(
        _mm_and_si128(mask_lower, cells),
        three
    ));
    uint32_t const higher16 = _mm_movemask_epi8(_mm_cmpeq_epi8(
        _mm_andnot_si128(mask_lower, cells),
        _mm_slli_epi16(three, 4)
    ));

    return lower16 | (higher16 << 16);
}

static Remainder calcRemainder(Field::FieldPimpl const &grid, int32_t const batchIndex) {
    auto const buffer = grid.getBuffer(Field::FieldPimpl::bufCur);
    auto const base = buffer + grid.bufferPaddingLength() + batchIndex -  1;
    auto const top  = *(base - grid.rowLength);
    auto const prev = *(base);
    auto const next = *(base + grid.rowLength);
    auto const at = [](uint32_t const value, int const offset) {
        return (value >> offset) & 1;
    };

    return {
        uint16_t(
            (at(top, 30) + at(prev, 30) + at(next, 30))
            + ((at(top, 31) + at(prev, 31) + at(next, 31)) << 8)
        ),
        (bool) at(prev, 31)
    };
}


static void threadUpdateGrid(Field::GridData& data) {
    Timer<> t{};
    auto& grid = *data.grid.get();
    int32_t const width_grid = static_cast<int32_t>(grid.width);
    int32_t const width_int = static_cast<int32_t>(grid.rowLength);


    int32_t const startBatch = data.startBatch;
    int32_t const endBatch = data.endBatch;
    int32_t const startIndex = startBatch * cellsBatchLength;
    int32_t const endIndex = endBatch * cellsBatchLength; //may be out of bounds padding required

    int32_t i = startBatch;
    auto previousRemainder = calcRemainder(grid, i);
    uint64_t newGenWindow = 0;

    auto const buffer = grid.getBuffer(Field::FieldPimpl::bufCur) + grid.bufferPaddingLength();
    auto const rowLen = grid.rowLength;

    if (i < endBatch + 1) {
        auto const newGen = newGenerationBatched(previousRemainder, buffer + i, rowLen, previousRemainder/*out param*/);

        newGenWindow = (uint64_t(newGen) << 31);

        ++i;
    }

    for (auto const j_count = 32; (i + j_count) < endBatch + 1;) {
        for (uint32_t j = 0; j < j_count; ++j, ++i) {
            auto const newGen = newGenerationBatched(previousRemainder, buffer + i, rowLen, previousRemainder/*out param*/);

            newGenWindow = (newGenWindow >> 32) | (uint64_t(newGen) << 31);
            grid.getCellsActual_int(i - 1, Field::FieldPimpl::bufNext) = uint32_t(newGenWindow);

            //_mm_stream_si32((int*)&grid.getCellsActual_int<Field::FieldPimpl::buffer>(i_batch - 1), (int)uint32_t(newGenWindow));
        }

        if (data.interrupt_flag.load()) return;
    }

    for (; i < endBatch + 1; ++i) {
        auto const newGen = newGenerationBatched(previousRemainder, buffer + i, rowLen, previousRemainder/*out param*/);

        newGenWindow = (newGenWindow >> 32) | (uint64_t(newGen) << 31);
        grid.getCellsActual_int(i - 1, Field::FieldPimpl::bufNext) = uint32_t(newGenWindow);
    }

    if (data.interrupt_flag.load()) return;


    if (grid.edgeCellsOptimization == false) {
        auto const setBufferCellAt = [&grid](uint32_t index, FieldCell cell) -> void { grid.setCellAt(index, cell, Field::FieldPimpl::bufNext); };

        auto const fullWidth = grid.rowLength * cellsBatchLength;
        uint32_t const height = grid.height;
        int32_t  const lastElement = width_grid - 1;

        const uint32_t startRow = startIndex / fullWidth;
        const uint32_t endRow = misc::min<uint32_t>(misc::intDivCeil(endIndex + 1, fullWidth), height);

        const auto isCell = [&grid](int32_t index) -> uint8_t {
            return grid.cellAt_grid(index);
        };
        const auto cellAt = [&grid](int32_t index) -> FieldCell {
            return grid.cellAt_grid(index);
        };

        {
            const int32_t row = startRow * width_grid;

            uint8_t //top/cur/bot + first/second/last
                tf = isCell(row + -width_grid),
                ts = isCell(row + -width_grid + 1),
                tl = isCell(row + -width_grid + lastElement),
                cf = isCell(row + 0),
                cs = isCell(row + 1),
                cl = isCell(row + lastElement);

            for (uint32_t rowIndex = startRow; rowIndex < endRow; rowIndex++) {
                const int32_t row = rowIndex * width_grid;

                uint8_t
                    bf = isCell(row + width_grid),
                    bs = isCell(row + width_grid + 1),
                    bl = isCell(row + width_grid + lastElement);
                //first element
                {
                    const int index = row;
                    const auto curCell = cellAt(index);
                    const unsigned char    topRowNeighbours = tf + ts + tl;
                    const unsigned char bottomRowNeighbours = bf + bs + bl;
                    const unsigned char aliveNeighbours = topRowNeighbours + cl + cs + bottomRowNeighbours;

                    setBufferCellAt(index, fieldCell::nextGeneration(curCell, aliveNeighbours));
                }

                tf = cf;
                ts = cs;
                tl = cl;
                cf = bf;
                cs = bs;
                cl = bl;
            }
        }
        if (data.interrupt_flag.load()) return;

        {
            const int32_t row = startRow * width_grid;

            uint8_t //top/cur/bot + first/last/pre-last
                tf = isCell(row + -width_grid),
                tl = isCell(row + -width_grid + lastElement),
                tp = isCell(row + -width_grid + lastElement - 1),
                cf = isCell(row + 0),
                cl = isCell(row + lastElement),
                cp = isCell(row + lastElement - 1);

            for (uint32_t rowIndex = startRow; rowIndex < endRow; rowIndex++) {
                const int32_t row = rowIndex * width_grid;

                uint8_t
                    bf = isCell(row + width_grid),
                    bp = isCell(row + width_grid + lastElement - 1),
                    bl = isCell(row + width_grid + lastElement);

                { //last element
                    const int index = row + lastElement;
                    const auto curCell = cellAt(index);
                    const uint8_t    topRowNeighbours = tp + tl + tf;
                    const uint8_t bottomRowNeighbours = bp + bl + bf;
                    const uint8_t     aliveNeighbours = topRowNeighbours + cp + cf + bottomRowNeighbours;

                    setBufferCellAt(index, fieldCell::nextGeneration(curCell, aliveNeighbours));
                }

                tf = cf;
                tp = cp;
                tl = cl;
                cf = bf;
                cp = bp;
                cl = bl;
            }
        }
        if (data.interrupt_flag.load()) return;
    }

    data.gridUpdate.add(t.elapsedTime());

    Timer<> t2{};
     
    data.buffer_output->write(FieldModification{ data.startBatch, data.endBatch - data.startBatch, &grid.getCellsActual_int(startBatch, Field::FieldPimpl::bufNext) });

    data.bufferSend.add(t2.elapsedTime());

    const uint32_t j = 1 << (data.task__iteration % (((grid.width-1) % 32) + 1));
    const uint32_t zero = 0;
    if(false) {
        auto output = data.buffer_output->batched();
        output->write(FieldModification{ 0, 1, &zero });
        output->write(FieldModification{ 1, 1, &j });
        output->write(FieldModification{ 2, 1, &zero });
    }

    data.generationUpdated();
}

uint32_t Field::width_actual() const {
    return gridPimpl->rowLength * cellsBatchLength;
}

uint32_t Field::size_bytes() const {
    return gridPimpl->gridLength() * cellsBatchSize;
}

uint32_t* Field::rawData() const {
    return &this->gridPimpl->getCellsActual_int(0);
}


Field::Field(
    const uint32_t gridWidth, const uint32_t gridHeight, const size_t numberOfTasks_,
    std::function<std::unique_ptr<FieldOutput>()> current_outputs, 
    std::function<std::unique_ptr<FieldOutput>()> buffer_outputs
) :
    gridPimpl{ new FieldPimpl(gridWidth, gridHeight) },
    isStopped{ false },
    current_output{ current_outputs() },
    buffer_output{ buffer_outputs() },
    numberOfTasks(numberOfTasks_),
    gridTasks{ new std::unique_ptr<Task<GridData>>[numberOfTasks_] },
    interrupt_flag{ false },
    indecesToBrokenCells{ }
{
    assert(numberOfTasks >= 1);

    //grid update == 50 ms, buffer send == 4 ms.In my case
    //4 / 50 ~= 1 / 20 
    const double fraction = .05; // fractionOfWorkTimeToSendData
    static_assert(true, ""/*
    find amountOfWorkT1. given fraction, workload (total batches)
    amountOfWorkT1 - percent of work;

    timeForTask1 = amountOfWorkT1 + amountOfWorkT1 * fraction;
    timeForTask2 = timeForTask1 + timeForTask1 * fraction;
    timeForTask3 = timeForTask2 + timeForTask2 * fraction;
    ...

                                                                v  total time
    timeForTask1 + timeForTask2 + timeForTask3 + ... = workload * (1.0 + fraction);
    amountOfWorkT1 * (1.0 + fraction)
        + (amountOfWorkT1 * (1.0 + fraction)) * (1.0  + fraction)
        + ...                                                       = workload * (1.0 + fraction);
    amountOfWorkT1 * (1.0 + fraction) + amountOfWorkT1 * (1.0  + fraction)^2 + ...  + amountOfWorkT1 * (1.0  + fraction)^numberOfTasks = workload * (1.0 + fraction);

    amountOfWorkT1 * ((1.0 + fraction) + (1.0  + fraction)^2 + ... + (1.0  + fraction)^numberOfTasks) = workload * (1.0 + fraction);

    amountOfWorkT1 = workload / ( ( (1.0 + fraction)^numberOfTasks - 1 ) / fraction );
    amountOfWorkT1 = workload /   ( (1.0 + fraction)^numberOfTasks - 1 ) * fraction  ;
    */);

    const auto createGridTask = [this, &buffer_outputs](const uint32_t index, const uint32_t startBatch, const uint32_t endBatch) -> void {
        glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);

        gridTasks.get()[index] = std::unique_ptr<Task<GridData>>(
            new Task<GridData>{
                threadUpdateGrid,
                
                index,
                this->gridPimpl,
                this->interrupt_flag,
                startBatch,
                endBatch,
                buffer_outputs() //getting output    
            }
        );
    };

    uint32_t batchesBefore = 0;
    uint32_t remainingBatches = gridPimpl->gridLength();

    const double amountOfWorkT1 = remainingBatches / (pow((1.0 + fraction), numberOfTasks) - 1) * fraction;

    double currentTaskWork = amountOfWorkT1;

    for (uint32_t i = 0; i < numberOfTasks - 1; i++) {
        const uint32_t numberOfBatches = misc::min(remainingBatches, static_cast<uint32_t>(ceil(currentTaskWork)));
        ::std::cout << numberOfBatches << std::endl;

        createGridTask(i, batchesBefore, batchesBefore + numberOfBatches);

        batchesBefore += numberOfBatches;
        remainingBatches -= numberOfBatches;

        currentTaskWork *= (1.0 + fraction);
    }

    {//last
        const uint32_t numberOfBatches = remainingBatches;
        ::std::cout << numberOfBatches << std::endl;

        createGridTask(numberOfTasks - 1, batchesBefore, batchesBefore + numberOfBatches);

        batchesBefore += numberOfBatches;
        remainingBatches -= numberOfBatches;
    }
}


Field::~Field() = default;

void Field::fill(const FieldCell cell) {
    interrupt_flag.store(true);
    waitForGridTasks();
    interrupt_flag.store(false);

    gridPimpl->fill(cell);

    indecesToBrokenCells.clear();

    current_output->write(FieldModification{ 0, static_cast<uint32_t>(gridPimpl->gridLength()), &gridPimpl->getCellsActual_int(0) });
    

    deployGridTasks();
}

FieldCell Field::cellAtIndex(const uint32_t index) const {
    return gridPimpl->cellAt_grid(normalizeIndex(index));
}

void Field::setCellAtIndex(const uint32_t index, FieldCell cell) {
    Cell const cell_{ cell, normalizeIndex(index) };
    setCells(&cell_, 1);
}

void Field::setCells(Cell const* const cells, size_t const count) {
    if (isStopped) {
        for (size_t i = 0; i < count; ++i) {
            auto const cell = cells[i];
            gridPimpl->setCellAt(normalizeIndex(cell.index), cell.cell);
        }
    }
    else {
        std::vector<uint32_t> cells_indeces_actual_int{};
        //std::vector<uint32_t> cells_int{};

        for (size_t i = 0; i < count; ++i) {
            auto const cell = cells[i];
            auto const index = normalizeIndex(cell.index);
            const auto index_actual_int{ gridPimpl->cellI2BatchI(index) };

            if (std::find(
                cells_indeces_actual_int.begin(), cells_indeces_actual_int.end(),
                index_actual_int
            ) == cells_indeces_actual_int.end()) {
                cells_indeces_actual_int.push_back(index_actual_int);
            }

            gridPimpl->setCellAt(index, cell.cell);
            indecesToBrokenCells.push_back(index);
        }

        for (uint32_t index_actual_int : cells_indeces_actual_int) {

            current_output->write(FieldModification{ index_actual_int, 1, &gridPimpl->getCellsActual_int(index_actual_int) });
        }
    }
}

static FieldCell updatedCell(const int32_t index, const std::unique_ptr<Field::FieldPimpl>& cellsGrid) {
    const auto width_grid = cellsGrid->width;
    const auto height = cellsGrid->height;
    uint8_t cell = cellsGrid->cellAt_grid(index);

    uint32_t aliveNeighbours = 0;

    for (int yo = -1; yo <= 1; yo++) {
        for (int xo = -1; xo <= 1; xo++) {
            if (xo != 0 || yo != 0) {
                int x = ((index + xo) + width_grid) % width_grid,
                    y = (((index / width_grid) + yo) + height) % height;
                int offsetedIndex = x + width_grid * y;
                if (cellsGrid->cellAt_grid(offsetedIndex))
                    aliveNeighbours++;
            }
        }
    }

    return fieldCell::nextGeneration(cell, aliveNeighbours);
}

bool Field::tryFinishGeneration() {
    if(!isStopped) for(uint32_t i = 0; i < numberOfTasks; i++) {
        if(!gridTasks.get()[i]->resultReady()) return false;
    }
    interrupt_flag.store(false);

    if (indecesToBrokenCells.size() > 0) {
        std::vector<uint32_t> repairedCells_actual_int{};

        for (auto it = indecesToBrokenCells.begin(); it != indecesToBrokenCells.end(); it++) {
            const uint32_t index = *it;
            const vec2i coord = this->indexAsCoord(index);

            for (int yo = -1; yo <= 1; yo++) {
                for (int xo = -1; xo <= 1; xo++) {
                    //int x = ((index + xo) + width()) % width(),
                    //    y = (((index / width()) + yo) + height()) % height();
                    auto offsetedIndex = this->coordAsIndex(coord + vec2i(xo, yo));//this->normalizeIndex(index + xo + gridPimpl->width_grid * yo);//x + width() * y;

                    const auto index_actual_int{ gridPimpl->cellI2BatchI(offsetedIndex) };

                    if (
                        std::find(
                            repairedCells_actual_int.begin(), repairedCells_actual_int.end(), 
                            index_actual_int
                        ) == repairedCells_actual_int.end()
                    ) repairedCells_actual_int.push_back(index_actual_int);
                }
            }
        }

        gridPimpl->fixField();

        std::unique_ptr<FieldOutput> output = buffer_output->batched();
        auto& field = *this->gridPimpl.get();
        auto const rowLen = field.rowLength;
        int32_t width_actual = field.rowLength * cellsBatchLength;
        int32_t width_grid = field.width;
        auto const buffer = field.getBuffer(Field::FieldPimpl::bufCur) + field.bufferPaddingLength();
        for (uint32_t index_actual_int : repairedCells_actual_int) {
            bool const isEdgeOpt = gridPimpl->edgeCellsOptimization;

            auto const colIndex = int(index_actual_int % rowLen);
            auto const specialFirstCell = !isEdgeOpt && (colIndex == 0);
            auto const specialLastCell = !isEdgeOpt && (colIndex == rowLen - 1);

            auto previousRemainder = ([&]() -> Remainder {
                    if(specialFirstCell) return {} /* we don't need to fill the remainder as it affects only the first two cells in the batch, one of them is from the previous row, and the second one is recalculated later*/;
                    else return calcRemainder(field, index_actual_int);
            })();

            uint64_t newGenerationWindow = newGenerationBatched(previousRemainder, buffer + index_actual_int, rowLen, previousRemainder);


            if(specialLastCell) {
                newGenerationWindow = (newGenerationWindow >> 1);
            }
            else {
                newGenerationWindow =
                    (newGenerationWindow >> 1) |
                    (uint64_t(newGenerationBatched(previousRemainder, buffer + index_actual_int + 1, rowLen, previousRemainder)) << (cellsBatchLength - 1));
            }

            uint32_t newGeneration = static_cast<uint32_t>(newGenerationWindow);

            auto const startRowCellIndex = index_actual_int / rowLen * width_grid;
            if(specialFirstCell) {
                newGeneration = (newGeneration & ~uint32_t(1)) | (updatedCell(startRowCellIndex, gridPimpl));
            }

            if(specialLastCell) {
                const auto lastCellCol = misc::mod(field.width-1, cellsBatchLength);
                newGeneration = (newGeneration & ~(uint32_t(1) << lastCellCol)) 
                    | (uint32_t(updatedCell(startRowCellIndex + field.width - 1, gridPimpl)) << lastCellCol);
            }


            gridPimpl->getCellsActual_int(index_actual_int, Field::FieldPimpl::bufNext) = newGeneration;
            output->write(FieldModification{ index_actual_int, 1, &newGeneration });
        }
        indecesToBrokenCells.clear();
    }

    return true;
}

void Field::startNewGeneration() {
    gridPimpl->swapBuffers();
    startCurGeneration();
}

void Field::startCurGeneration() {
    gridPimpl->fixField();
    isStopped = false;
    deployGridTasks();
}

void Field::waitForGridTasks() {
    if (isStopped) return;
    for (uint32_t i = 0; i < numberOfTasks; i++) {
        gridTasks.get()[i]->waitForResult();
    }
}
void Field::deployGridTasks() {
    if(isStopped) {
        std::cerr << "trying to start task when `isStopped` is set\n";
        return;
    }
    for(uint32_t i = 0; i < numberOfTasks; i++) {
        gridTasks.get()[i]->start();
    }
}

uint32_t Field::width() const {
    return gridPimpl->width;
}
uint32_t Field::height() const {
    return gridPimpl->height;
}
uint32_t Field::size() const {
    return gridPimpl->width * gridPimpl->height;
}
