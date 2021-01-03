#include <iostream>
#include <bits/stdc++.h>

using namespace std;
using namespace chrono;

// Integer types.
using Board = long long;
using Zones = long long;
using Move = uint;

// Board sizes.
const int W = 5;
const int L = W + 1;
const int N = 2 * W * L;
const int M = W + L;

// Constant boards.
const Board ONE = 1;
const Board EMPTY = 0;
const Board FULL = 1152921504606846975;
const Board LEFT = 563224965611552;
const Board NOT_LEFT = LEFT ^FULL;
const Board RIGHT = 18023198899569664;
const Board NOT_RIGHT = RIGHT ^FULL;
const Board TOP = 31;
const Board BOTTOM = 1116892707587883008;
const Board VERTICAL = 35483172833527776;
const Board HORIZONTAL = 1117438331773319199;
const Board BORDER = (LEFT & VERTICAL) | (RIGHT & VERTICAL) | TOP | BOTTOM;

// Constant io.
const string START = "Start";
const string QUIT = "Quit";
const string HLINE = "▁▁";
const string HDASH = "__";
const string VLINE = "┃";
const string VDASH = "┊";


// Combined types.
struct State {
    Board b;
    Zones z;

    ~State()= default;

    // Required for hash map.
    bool operator==(const State &s) const {
        return b == s.b && z == s.z;
    }
};

namespace std {

    template<>
    struct [[maybe_unused]] hash<State> {
        std::size_t operator()(const State &s) const {
            return hash<Board>()(s.b) ^ hash<Zones>()(s.z);
        }
    };
}

// MCTS state.
using Visits = unordered_map<State, ulong>;
using Plays = unordered_map<State, unordered_map<Move, ulong>>;
using Wins = unordered_map<State, unordered_map<Move, ulong>>;
using Children = unordered_map<State, unordered_map<Move, State>>;


class MctsState {
public:
    float winChances(State s, Move m) {
        return ((float) getWins(s, m)) / ((float) plays[s][m]);
    }

    void updateLeaf(State s) {
        visits[s]++;
    }

    void update(State s, Move m, uint win) {
        visits[s]++;
        plays[s][m]++;
        wins[s][m] += win;
    }

    ulong getVisits(State s) {
        return visits[s];
    }

    ulong getPlays(State s, Move m) {
        auto got = plays.find(s);
        if (got == plays.end()) {
            return 0;
        } else {
            return got->second[m];
        }
    }

    ulong getWins(State s, Move m) {
        auto got = wins.find(s);
        if (got == wins.end()) {
            return 0;
        } else {
            return got->second[m];
        }
    }

    void setChildren(State s, unordered_map<Move, State> &moves) {
        children[s] = moves;
    }

    bool hasChildren(State s) {
        return !(children.find(s) == children.end());
    }

    unordered_map<Move, State> getChildren(State s) {
        auto got = children.find(s);
        if (got == children.end()) {
            return unordered_map<Move, State>{};
        } else {
            return got->second;
        }
    }

private:
    Visits visits = Visits(10000);
    Plays plays = Plays(10000);
    Wins wins = Wins(10000);
    Children children = Children(10000);
};

// Conversion tables.
string INT_TO_POSITION[N] = {"A1h", "A2h", "A3h", "A4h", "A5h",
                             "A1v", "A2v", "A3v", "A4v", "A5v", "A6v",
                             "B1h", "B2h", "B3h", "B4h", "B5h",
                             "B1v", "B2v", "B3v", "B4v", "B5v", "B6v",
                             "C1h", "C2h", "C3h", "C4h", "C5h",
                             "C1v", "C2v", "C3v", "C4v", "C5v", "C6v",
                             "D1h", "D2h", "D3h", "D4h", "D5h",
                             "D1v", "D2v", "D3v", "D4v", "D5v", "D6v",
                             "E1h", "E2h", "E3h", "E4h", "E5h",
                             "E1v", "E2v", "E3v", "E4v", "E5v", "E6v",
                             "F1h", "F2h", "F3h", "F4h", "F5h"};

// Board functions.
Board put(Board b, Move m) {
    if ((b >> m & 1) > 0) {
        throw invalid_argument("Cannot place wall at " + to_string(m)
                               + " because it already has a wall");
    }
    return b | ONE << m;
}

Board proceed(Board b) {
    return b
           | ((b & VERTICAL & NOT_RIGHT) << 1)
           | ((b & VERTICAL & NOT_LEFT) >> 1)
           | ((b & NOT_LEFT) << W)
           | ((b & NOT_RIGHT) >> W)
           | ((b & NOT_RIGHT) << (W + 1))
           | ((b & NOT_LEFT) >> (W + 1))
           | ((b & HORIZONTAL) << M)
           | ((b & HORIZONTAL) >> M);
}

Board fillEnclosed(Board b) {
    Board bb = (b ^ FULL) & BORDER;
    Board previous;
    do {
        previous = bb;
        bb = proceed(bb) & (b ^ FULL);
    } while (bb != previous);

    return b | (bb ^ FULL);
}

// Zone functions.
uint countOnes(Board b) {
    return ((bitset<64>) b).count();
}

uint countEnclosedTiles(Board b) {
    b = fillEnclosed(b);
    Board shiftBack = b & (b >> M) & (b >> W) & (b >> (W + 1));
    return countOnes(b & shiftBack & HORIZONTAL & ~BOTTOM);
}

uint countAddedZones(Board b, Board bb) {
    return countEnclosedTiles(bb) - countEnclosedTiles(b);
}

list<uint> listZones(Zones z) {
    list<uint> zoneSizes;
    auto zBit = (bitset<64>) z;
    for (ulong i = 0; i < zBit.size(); i++) {
        if (zBit[i]) {
            zoneSizes.push_back(i);
        }
    }
    return zoneSizes;
}

Zones updateZones(Zones z, uint zone) {
    if (zone == 0) {
        return z;
    } else {
        return z | 1 << zone;
    }
}

// State functions.

State play(State s, Move m) {
    Board bb = fillEnclosed(put(s.b, m));
    uint addedZones = countAddedZones(s.b, bb);
    State ss {bb, updateZones(s.z, addedZones)};
    return ss;
}

unordered_map<Move, State> validMoves(State s) {
    unordered_map<Move, State> moves;
    bitset<N> bits = s.b;
    for (int i = 0; i < N; ++i) {
        if (bits[i] == 0) {
            Board bb = fillEnclosed(put(s.b, i));
            uint addedZones = countAddedZones(s.b, bb);

            if (((1 << addedZones) & s.z) == 0) {
                State ss {bb, (Zones) updateZones(s.z, addedZones)};
                moves[i] = ss;
            }
        }
    }
    return moves;
}

tuple<bool, Move, State> randomMove(State s) {
    bitset<N> bBits = s.b;
    uint shift = random() % N;
    for (uint i = 0; i < N; i++) {
        uint j = (i + shift) % N;
        if (bBits[j] == 0) {
            Board bb = fillEnclosed(put(s.b, j));
            uint addedZones = countAddedZones(s.b, bb);

            if (((1 << addedZones) & s.z) == 0) {
                State ss {bb, (Zones) updateZones(s.z, addedZones)};
                return tuple<bool, Move, State>(true, j, ss);
            }
        }
    }
    return tuple<bool, Move, State>(false, 0, State {});
}

// mcts functions.

ulong simulate(State s) {
    auto valid_move_state = randomMove(s);
    if (get<0>(valid_move_state)) {
        return 1 - simulate(get<2>(valid_move_state));
    } else {
        return 0;
    }
}

float score(ulong w, ulong n, ulong nParent, float c, float eps = 0.0001) {
    float n_eps = n + eps;
    float n_parent_eps = nParent + 1.0F;

    return ((float) w / n_eps)
           + c * sqrtf(std::log(n_parent_eps) / n_eps);
}

pair<Move, State> bestMove(MctsState &mcts_state, State s, float c) {
    unordered_map<Move, State> children = mcts_state.getChildren(s);
    pair<Move, State> best_move_state;
    float best_score = -INFINITY;
    for (auto move_state: children) {
        const Move m = move_state.first;

        ulong wins = mcts_state.getWins(s, m);
        ulong plays = mcts_state.getPlays(s, m);
        ulong visits = mcts_state.getVisits(s);
        const float m_score = score(wins,
                                    plays,
                                    visits, c);
        if (m_score > best_score) {
            best_move_state = move_state;
            best_score = m_score;
        }
    }
    return best_move_state;
}


uint search(MctsState &mcts_state, State s, float c) {
    if (mcts_state.hasChildren(s)) {
        pair<Move, State> move_and_state = bestMove(mcts_state, s, c);
        uint win = 1 - search(mcts_state, move_and_state.second, c);

        mcts_state.update(s, move_and_state.first, win);

        return win;
    } else {
        // s not in mcts_state.children
        unordered_map<Move, State> moves = validMoves(s);
        if (moves.empty()) {
            mcts_state.updateLeaf(s);
            return 0;
        } else {
            mcts_state.setChildren(s, moves);

            // simulate.
            pair<Move, State> move_and_state = bestMove(mcts_state, s, c);
            uint win = 1 - simulate(move_and_state.second);

            // update leaf node.
            mcts_state.update(s, move_and_state.first, win);

            return win;
        }
    }
}


pair<Move, State>
mcts(MctsState &mcts_state, State s, float c, float duration_secs) {
    auto start = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(
            high_resolution_clock::now() - start);

    do {
        search(mcts_state, s, c);
        duration = duration_cast<microseconds>(
                high_resolution_clock::now() - start);
    } while ((duration.count() / 1000.0 / 1000.0) < duration_secs);

    auto children = mcts_state.getChildren(s);
    pair<Move, State> best_move_state;
    uint max_plays = 0;
    for (auto move_state : children) {
        ulong plays = mcts_state.getPlays(s, move_state.first);
        if (plays > max_plays) {
            max_plays = plays;
            best_move_state = move_state;
        }
    }
    return best_move_state;
}


// Game functions.

string read() {
    string input;
    std::cin >> input;
    return input;
}

void write(const string &s) {
    std::cout << s << std::endl;
    std::cout.flush();
}

void log(const string &s) {
    std::cerr << s;
}

void logMcts(MctsState mcts_state, State s, Move m) {
    float win_chances = mcts_state.winChances(s, m);
    log("Chances of winning are " + to_string(win_chances) + " after "
        + to_string(mcts_state.getVisits(s)) + " moves\n");
}

void logBoard(Board b) {
    auto bBit = (bitset<N>) b;
    log("\nBoard " + to_string(b) + "\nBin " + bBit.to_string() + "\n");
    for (int i = 0; i < W; i++) {
        if (bBit[i]) {
            log(" " + HLINE);
        } else {
            log(" " + HDASH);
        }
    }
    log(" \n");
    for (int i = 0; i < W; i++) {
        for (int j = 0; j < W; j++) {
            if (bBit[W + M * i + j]) {
                log(VLINE);
            } else {
                log(VDASH);
            }
            if (bBit[W + M * i + j + L]) {
                log(HLINE);
            } else {
                log(HDASH);
            }
        }
        if (bBit[W + M * i + W]) {
            log(VLINE);
        } else {
            log(VDASH);
        }
        log("\n");
    }
}

void logState(State s, uint verbose = 2) {
    if (verbose >= 2) {
        logBoard(s.b);
    }
    if (verbose >= 1) {
        list<uint> zoneSizes = listZones(s.z);
        log("Enclosed zones: ");
        for (uint i: zoneSizes) {
            log(to_string(i) + ", ");
        }
        log("\n");
    }
    if (verbose >= 0) {
        unordered_map<Move, State> moves = validMoves(s);
        log("Possible moves: ");
        for (pair<Move, State> moveState: moves) {
            log(INT_TO_POSITION[moveState.first] + ", ");
        }
        log("\n");
    }
}

float timeBudget(State s, float total_time) {
    float time_left = 29.0F - total_time;
    float moves_left = validMoves(s).size();
    float max_move_time = time_left / 2.0F;
    float min_move_time = time_left / (moves_left + 1);

    if (moves_left < 15) {
        return min(min_move_time * 4.0F, max_move_time);
    } else if (moves_left < 30) {
        return min(min_move_time * 2.0F, max_move_time);
    } else {
        return min(min_move_time * 0.1F, max_move_time);
    }
}

uint positionToInt(const string& pos) {
    for (ulong i = 0; i < sizeof(INT_TO_POSITION); i++) {
        if (INT_TO_POSITION[i] == pos) {
            return i;
        }
    }
    throw invalid_argument("Cannot convert '" + pos + "' to integer");
}

void game() {
    float total_ms = 0;
    State s{};

    string text = read();
    while (text != QUIT) {
        auto start_time = high_resolution_clock::now();

        if (text != START) {
            s = play(s, positionToInt(text));
        }
        float duration = timeBudget(s, total_ms / 1000.0F);
        log("Spending " + to_string(duration) + "s of "
            + to_string(30.0F - total_ms / 1000.0F) + "s that are left.\n");

        MctsState mcts_state = {};
        auto move_state = mcts(mcts_state, s, sqrtf(2.0), duration);
        logMcts(mcts_state, s, move_state.first);
//        mcts_state.clean(s);
        s = move_state.second;

        logState(s, 1);
        write(INT_TO_POSITION[move_state.first]);

        auto end_time = high_resolution_clock::now();
        total_ms += duration_cast<milliseconds>(end_time - start_time).count();
        text = read();
    }
}

// Testing.
const Board ENCLOSED_SQUARE = 52882452480;
const Board SINGLE_TILE = 17993564160;

// Testing board.
bool test_put_peg() {
    Board b = put(EMPTY, 0);
    return b == 1;
}

bool test_putting_peg_twice_does_not_change_board() {
    Board b = put(EMPTY, 10);
    try {
        put(b, 10);
    } catch(...) {
        return true;
    }
    return false;
}

bool test_putting_all_pegs() {
    Board b = EMPTY;
    for (Move m = 0; m < N; m++) {
        b = put(b, m);
    }
    return b == FULL;
}

bool testFillEnclosed() {
    const Board b = fillEnclosed(ENCLOSED_SQUARE);
    return b > ENCLOSED_SQUARE;
}

bool testFillEmptyIsEmpty() {
    const Board b = fillEnclosed(EMPTY);
    return b == EMPTY;
}

bool testFiveRandomPegsEncloseNoOtherPegs() {
    Board b = EMPTY;
    for (int i = 0; i < 5; i++) {
        const Move m = random() % N;
        b = put(b, m);
    }
    Board bb = fillEnclosed(b);
    return b == bb;
}

bool testFullIsFull() {
    Board b = fillEnclosed(FULL);
    return b == FULL;
}

bool testFillBoarderIsFull() {
    Board b = fillEnclosed(BORDER);
    return b == FULL;
}

// Testing zones.
bool testEmptyBoardHasZeroOnes() {
    return countEnclosedTiles(EMPTY) == 0;
}

bool testFullBoardEnclosesAll() {
    return countEnclosedTiles(FULL) == W * W;
}

bool testBorderEnclosesAll() {
    return countEnclosedTiles(BORDER) == W * W;
}

bool testEnclosedSquareHas4Tiles() {
    return countEnclosedTiles(ENCLOSED_SQUARE) == 4;
}

bool testEnclosedTileHas1Tile() {
    return countEnclosedTiles(SINGLE_TILE) == 1;
}

bool testAddedZoneSize() {
    Board b = 869832927612085311;
    Board bb = fillEnclosed(put(b, 9));  // A5v
    return countAddedZones(b, bb) == 1;
}

// Test state.

// Regression test for https://www.codecup.nl/showgame.php?ga=164794
bool testInvalidMove() {
    Board b = 1152894691084558198;
    auto s = State {b, 30};
    State s1 = play(s, 20);
    State s2 = play(s, 21);
    State s3 = play(s, 32);
    State s4 = play(s, 37);
    State s5 = play(s, 43);

    unordered_map<Move, State> moves = validMoves(s);
    return moves.size() == 5
           && moves[20] == s1
           && moves[21] == s2
           && moves[32] == s3
           && moves[37] == s4
           && moves[43] == s5;
}

// Test mcts.

bool testSimulate() {
    uint value = simulate(State {});
    return value == 0 || value == 1;
}

bool testSpeedSimulate() {
    auto start = high_resolution_clock::now();
    for (int i = 0; i < 1000; i++) {
        simulate(State {});
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
//    log(to_string(duration.count()));
    return duration.count() < 50000;  // 50ms, 50 microseconds per simulation
}

bool testScore() {
    return score(10, 20, 400, 2) < score(11, 20, 400, 2)
           && score(10, 20, 400, 2) < score(10, 19, 400, 2)
           && score(10, 20, 400, 1) < score(10, 20, 400, 2);
}

bool testBestMoveWithThreeSquares() {
    State s{1152921504602650622, 0};
    MctsState mcts_state;
    auto move_state = mcts(mcts_state, s, 2.0, 0.05);

    return move_state.first == 0
           && mcts_state.winChances(s, move_state.first) > 0.9;
}

bool testBestMoveWithTwoSquares() {
    State s{1152921504606844926, 0};
    MctsState mcts_state;

    auto move_state = mcts(mcts_state, s, 2.0, 0.05);

    // both moves are a win
    return mcts_state.winChances(s, move_state.first) > 0.9;
}

bool testLosingWithThreeSquares() {
    State s{1152921504602650622, 1 << 3};
    MctsState mcts_state;

    auto move_state = mcts(mcts_state, s, 2.0, 0.05);

    return mcts_state.winChances(s, move_state.first) < 0.1;
}

bool testWinAgainstRandomized() {
    State s{0, 0};
    MctsState mcts_state;
    bool is_mcts = true;

    unordered_map<Move, State> moves;
    do {
        Move m;
        if (is_mcts) {
            auto move_state = mcts(mcts_state, s, 2.0, 0.05);
            m = move_state.first;
        } else {
            auto valid_move_state = randomMove(s);
            m = get<1>(valid_move_state);
        }
        s = play(s, m);
        moves = validMoves(s);
        is_mcts = !is_mcts;
    } while (!moves.empty());

    return !is_mcts;
}

[[maybe_unused]] void test() {
    std::cout << "Tests" << std::endl;
    std::cout << test_put_peg()
              << " test_put_peg" << std::endl;
    std::cout << test_putting_peg_twice_does_not_change_board()
              << " test_putting_peg_twice_does_not_change_board" << std::endl;
    std::cout << test_putting_all_pegs()
              << " test_putting_all_pegs" << std::endl;
    std::cout << testFillEnclosed()
              << " testFillEnclosed()" << std::endl;
    std::cout << testFillEmptyIsEmpty()
              << " testFillEmptyIsEmpty" << std::endl;
    std::cout << testFiveRandomPegsEncloseNoOtherPegs()
              << " testFiveRandomPegsEncloseNoOtherPegs" << std::endl;
    std::cout << testFullIsFull()
              << " testFullIsFull" << std::endl;
    std::cout << testFillBoarderIsFull()
              << " testFillBoarderIsFull" << std::endl;
    std::cout << testEmptyBoardHasZeroOnes()
              << " testEmptyBoardHasZeroOnes" << std::endl;
    std::cout << testFullBoardEnclosesAll()
              << " testFullBoardEnclosesAll" << std::endl;
    std::cout << testBorderEnclosesAll()
              << " testBorderEnclosesAll" << std::endl;
    std::cout << testEnclosedSquareHas4Tiles()
              << " testEnclosedSquareHas4Tiles" << std::endl;
    std::cout << testEnclosedTileHas1Tile()
              << " testEnclosedTileHas1Tile" << std::endl;
    std::cout << testAddedZoneSize()
              << " testAddedZoneSize" << std::endl;
    std::cout << testInvalidMove()
              << " testInvalidMove" << std::endl;
    std::cout << testSimulate()
              << " testSimulate" << std::endl;
    std::cout << testSpeedSimulate()
              << " testSpeedSimulate" << std::endl;
    std::cout << testScore()
              << " testScore" << std::endl;
    std::cout << testBestMoveWithThreeSquares()
              << " testBestMoveWithThreeSquares" << std::endl;
    std::cout << testBestMoveWithTwoSquares()
              << " testBestMoveWithTwoSquares" << std::endl;
    std::cout << testLosingWithThreeSquares()
              << " testLosingWithThreeSquares" << std::endl;
    std::cout << testWinAgainstRandomized()
              << " testWinAgainstRandomized" << std::endl;

}

int main() {
//    test();
    game();

    // todo move code into classes
    // todo reduce memory
    // todo figure out good simulate() for start of game

    return 0;
}
