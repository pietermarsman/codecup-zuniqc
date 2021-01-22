#include <iostream>
#include <bits/stdc++.h>

using namespace std;
using namespace chrono;

// Integer types.
using Board = long long;
using Zones = long;
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
class State {
public:
    Board b = 0;
    Zones z = 0;
    uint depth = 0;
    uint extra = 0;
    Board excluded = 0;

    // State functions.
    [[nodiscard]] optional<State> play(Move m) {
        if (m >= N) {
            return nullopt;
        }
        const auto maybe_b = put(b, m);
        if (!maybe_b.has_value()) {
            return nullopt;
        }

        Board bb = fillEnclosed(maybe_b.value());
        uint zone = countAddedZones(b, bb);
        if ((z >> zone & 1) == 1) {
            // zone already taken: check if move can always be excluded
            int smaller_zones = (1 << (zone + 1)) - 1;
            bool exclude = (smaller_zones & z) == smaller_zones;
            if (exclude) {
                excluded = excluded | (ONE << m);
            }

            return nullopt;
        }
        Zones zz = updateZones(z, zone);

        uint addedExtra = 0;
        if (zone > 1) {
            Board inner = bb & ~maybe_b.value();
            addedExtra = countOnes(inner) + 1 - zone;
        }

        return State{bb, zz, depth + 1, extra + addedExtra,
                     excluded};
    }

    // List all valid moves.
    [[nodiscard]] vector<pair<Move, State>> validMoves() {
        vector<pair<Move, State>> moves;
        bitset<N> bits = ~(b | excluded);
        for (long unsigned int i = bits._Find_first();
             i < bits.size(); i = bits._Find_next(i)) {
            auto maybe_s = play(i);
            if (maybe_s.has_value()) {
                moves.emplace_back(i, maybe_s.value());
            }
        }
        return moves;
    }

    // List single next move that is equal or larger to `nextMoveIndex`.
    [[nodiscard]] optional<pair<Move, State>> nextMove(ulong nextMoveIndex) {
        bitset<N> bits = ~(b | excluded);
        optional<pair<Move, State>> move_state = nullopt;
        if (nextMoveIndex == 0) {
            nextMoveIndex = bits._Find_first();
        }
        do {
            auto maybe_s = play(nextMoveIndex);
            if (maybe_s.has_value()) {
                 move_state = {nextMoveIndex, maybe_s.value()};
            }
            // always makes nextMoveIndex larger, unless nextMoveIndex is larger
            // than bits.size()
            nextMoveIndex = bits._Find_next(nextMoveIndex);
        } while (nextMoveIndex < bits.size() && !move_state.has_value());
        return move_state;
    }

    // Return single possible random move or `nullopt`.
    [[nodiscard]] optional<pair<Move, State>> randomMove() {
        bitset<N> bits = ~(b | excluded);
        uint shift = random() % N;
        // from `shift` to end
        for (long unsigned int i = bits._Find_next(shift);
             i < bits.size(); i = bits._Find_next(i)) {
            auto maybe_s = play(i);
            if (maybe_s.has_value()) {
                return pair<Move, State>{i, maybe_s.value()};
            }
        }
        // from start to `shift`
        for (uint i = bits._Find_first();
             i < shift; i = bits._Find_next(i)) {
            auto maybe_s = play(i);
            if (maybe_s.has_value()) {
                return pair<Move, State>{i, maybe_s.value()};
            }
        }
        return nullopt;
    }

    // List all zones by their size
    [[nodiscard]] vector<uint> listZones() const {
        vector<uint> zoneSizes;
        auto zBit = (bitset<64>) z;
        for (ulong i = 0; i < zBit.size(); i++) {
            if (zBit[i]) {
                zoneSizes.push_back(i);
            }
        }
        return zoneSizes;
    }

    // Required for hash map.
    bool operator==(const State &s) const {
        return b == s.b && z == s.z;
    }

    // Simulate a game by playing random moves.
    [[nodiscard]] ulong simulate() {
        auto maybe_state = randomMove();
        if (maybe_state.has_value()) {
            return 1 - maybe_state.value().second.simulate();
        } else {
            return 0;
        }
    }

private:
    // Board functions.
    static optional<Board> put(Board b, Move m) {
        if ((b >> m & 1) > 0) {
            return nullopt;
        }
        return b | ONE << m;
    }

    // Proceed board one step by enabling all walls one step away from an
    // already enabled wall.
    static Board proceed(Board b) {
        const Board bh = b & HORIZONTAL;
        const Board bnr = b & NOT_RIGHT;
        const Board bnl = b & NOT_LEFT;
        return b
               | ((bnr & VERTICAL) << 1)
               | ((bnl & VERTICAL) >> 1)
               | (bnl << W)
               | (bnr >> W)
               | (bnr << (W + 1))
               | (bnl >> (W + 1))
               | (bh << M)
               | (bh >> M);
    }

    // Enable all enclosed walls.
    static Board fillEnclosed(Board b) {
        Board bb = (b ^ FULL) & BORDER;
        Board previous;
        do {
            previous = bb;
            bb = proceed(bb) & (b ^ FULL);
        } while (bb != previous);

        return b | (bb ^ FULL);
    }

    // Count ones on board.
    static uint countOnes(Board b) {
        return ((bitset<64>) b).count();
    }

    // Count number of tiles that are enabled on all four edges.
    static uint countEnclosedTiles(Board b) {
        Board shiftBack = b & (b >> M) & (b >> W) & (b >> (W + 1));
        return countOnes(shiftBack & HORIZONTAL & ~BOTTOM);
    }

    // Count the number of added tiles with all four edges enabled.
    static uint countAddedZones(Board b, Board bb) {
        return countEnclosedTiles(bb) - countEnclosedTiles(b);
    }

    // Update the zones of the state.
    static Zones updateZones(Zones z, uint zone) {
        if (zone == 0) {
            return z;
        } else {
            return z | 1 << zone;
        }
    }
};


class Node {
public:
    Node *p = nullptr;
    Move m = 0;  // most recent move

    bool expanded = false;
    ulong nextMoveIndex = 0;
    vector<shared_ptr<Node>> children;

    ulong n = 0;
    float nLog = std::log(0);
    ulong w = 0;

    // Create root node of t; no p
    explicit Node() = default;

    // Create intermediate node of t
    explicit Node(Node *parent, Move move) {
        p = parent;
        m = move;
    }

    // Return best move according to the mcts scoring function.
    shared_ptr<Node> bestMove(float c) {
        shared_ptr<Node> best;
        float best_score = -INFINITY;
        for (const auto &child: children) {
            const float link_score = child->score(c);
            if (link_score > best_score) {
                best = child;
                best_score = link_score;
            }
        }
        return best;
    }

    // Compute score of node from perspective of parent.
    [[nodiscard]] float score(float c, float eps = 0.0001) const {
        float n_eps = n + eps;
        float wins = (float) n - w;

        return ((float) wins / n_eps) + c * std::sqrt(p->nLog / n_eps);
    }

    // Compute wins of node from perspective of parent.
    [[nodiscard]] ulong wins() const {
        return n - w;
    }

    // Compute win chances from perspective of parent.
    [[nodiscard]] float winChances() const {
        return (float) wins() / (float) n;
    }

    // Expand one option of node and return the new node, or `nullptr` if
    // it cannot be expanded further.
    pair<shared_ptr<Node>, State> expand(State s) {
        auto move_state = s.nextMove(nextMoveIndex);
        if (move_state.has_value()) {
            nextMoveIndex = move_state->first + 1;
            auto node = make_shared<Node>(this, move_state->first);
            children.push_back(node);
            return {node, move_state->second};
        } else {
            expanded = true;
            return {nullptr, State{}};
        }
    }

    // Back propagate win through all ancestors.
    void backPropagate(ulong win) {
        n++;
        nLog = std::log(n);
        w += win;
        if (p) {
            p->backPropagate(1 - win);
        }
    }
};


class Mcts {
public:
    explicit Mcts(State s, float constant, ulong expandAfterVisits = 1) {
        root = s;
        t = make_shared<Node>();
        c = constant;
        minExpand = expandAfterVisits;
    }

    // Execute mcts searches for fixed duration or number of iterations.
    shared_ptr<Node>
    mcts(float duration_secs, float max_iterations = INFINITY) {
        auto start = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(
                high_resolution_clock::now() - start);

        ulong i = 0;
        do {
            i++;
            search();
            duration = duration_cast<microseconds>(
                    high_resolution_clock::now() - start);
        } while ((duration.count() / 1000.0 / 1000.0) < duration_secs
                 && i < max_iterations);

        return getBestMove();
    }

    // Change the root of the tree to the state that belongs to the `move`.
    State updateRoot(Move move) {
        while (!t->expanded) {
            auto ret = t->expand(root);
        }
        Move m = 0;
        for (const auto &child: t->children) {
            if (child->m == move) {
                t = child;
                m = move;
                break;
            }
        }
        root = root.play(m).value();
        t->p = nullptr;
        return root;
    }

    // Compute time to spend on computing the next best move.
    [[nodiscard]] float timeBudget(float total_time) const {
        float time_left = 28.0F - total_time;
        float moves_left = 41.0F - root.depth;
        float max_move_time = time_left / 2.0F;
        float min_move_time = time_left / (moves_left + 1);

        if (moves_left < 15) {
            return min(min_move_time * 4.0F, max_move_time);
        } else if (moves_left < 27) {
            return min(min_move_time * 2.0F, max_move_time);
        } else {
            return min(min_move_time * 0.2F, max_move_time);
        }
    }

    // Get the tree.
    shared_ptr<Node> getTree() {
        return t;
    }

private:
    State root;
    shared_ptr<Node> t;
    float c;
    ulong minExpand;

    // Do single search, expanding node and simulate random playout.
    void search() const {
        shared_ptr<Node> n = t;
        State s = root;
        while (n->expanded && !n->children.empty()) {
            n = n->bestMove(c);
            s = s.play(n->m).value();
        }
        if (!n->expanded && n->n >= minExpand) {
            // stopped because it is not yet expanded, not because it has no
            // children.
            auto nn_ss = n->expand(s);
            if (nn_ss.first) {
                n = nn_ss.first;
                s = nn_ss.second;
            }
            // ignore else, because sometimes n does not have children
        }
        ulong win = s.simulate();
        n->backPropagate(win);
    }

    // Best move to play is the one with the most simulations
    shared_ptr<Node> getBestMove() {
        shared_ptr<Node> best_node;
        uint max_n = 0;
        for (const auto &child : t->children) {
            if (child->n > max_n) {
                max_n = child->n;
                best_node = child;
            }
        }
        return best_node;
    }
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

void logMcts(const shared_ptr<Node> &node) {
    log("Chances of winning are " + to_string(node->wins()) + " / "
        + to_string(node->n) + " = " + to_string(node->winChances()) + "\n");
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
        auto zoneSizes = s.listZones();
        log("Enclosed zones: ");
        for (uint i: zoneSizes) {
            log(to_string(i) + ", ");
        }
        log("\n");
        log("Extra enclosed walls: " + to_string(s.extra) + "\n");
    }
    if (verbose >= 0) {
        auto moves = s.validMoves();
        log("Possible moves: ");
        for (pair<Move, State> moveState: moves) {
            log(INT_TO_POSITION[moveState.first] + ", ");
        }
        log("\n");
    }
}

uint positionToInt(const string &pos) {
    for (ulong i = 0; i < sizeof(INT_TO_POSITION); i++) {
        if (INT_TO_POSITION[i] == pos) {
            return i;
        }
    }
    throw invalid_argument("Cannot convert '" + pos + "' to integer");
}

bool sureWin(const shared_ptr<Node> &node) {
    if ((node->n > 10000) && (node->winChances() > 0.8)) {
        return true;
    }
    return false;
}

void game(float c, ulong expandAfter) {
    float total_ms = 0;
    Mcts mcts = Mcts(State{}, c, expandAfter);
    bool i_am_winning = false;

    string text = read();
    while (text != QUIT) {
        auto start_time = high_resolution_clock::now();

        if (text != START) {
            // assuming input is always correct
            Move move = positionToInt(text);
            mcts.updateRoot(move);
        }
        float duration = mcts.timeBudget(total_ms / 1000.0F);
        log("Spending " + to_string(duration) + "s of "
            + to_string(30.0F - total_ms / 1000.0F) + "s that are left.\n");

        auto best_node = mcts.mcts(duration);
        State root = mcts.updateRoot(best_node->m);
        logMcts(mcts.getTree());
        logState(root, 1);

        if (!i_am_winning && sureWin(best_node)) {
            write(INT_TO_POSITION[best_node->m] + "!");
            i_am_winning = true;
        } else {
            write(INT_TO_POSITION[best_node->m]);
        }

        auto end_time = high_resolution_clock::now();
        total_ms += duration_cast<milliseconds>(end_time - start_time).count();
        text = read();
    }
}

// Testing board.
bool testPlayWall() {
    State s{EMPTY, 0};
    auto maybe_state = s.play(0);
    return maybe_state.value().b == 1;
}

bool testCannotPutWallTwiceOnSamePosition() {
    State s{EMPTY, 0};
    auto maybe_state = s.play(10);
    auto maybe_state2 = maybe_state.value().play(10);
    return !maybe_state2.has_value();
}

bool testPlaceAllBorderWalls() {
    optional<State> maybe_state{{EMPTY, 0}};
    // play border
    for (uint i = 0; i < W; i++) {
        maybe_state = maybe_state.value().play(i);
        maybe_state = maybe_state.value().play(N - i - 1);
        maybe_state = maybe_state.value().play(i * M + M - 1);
        maybe_state = maybe_state.value().play(i * M + W);
    }
    return maybe_state.value().b == FULL
           && maybe_state.value().z == 1 << W * W;
}

bool testPlayingOneTile() {
    State s{EMPTY, 0};
    s = s.play(0).value();
    s = s.play(W).value();
    s = s.play(W + 1).value();
    s = s.play(M).value();
    return (s.z = 1);
}

bool testPlayingSquareTile() {
    State s{EMPTY, 0};
    s = s.play(0).value();
    s = s.play(1).value();
    s = s.play(W).value();
    s = s.play(M + W).value();
    s = s.play(W + 2).value();
    s = s.play(M + W + 1).value();
    s = s.play(M + 1).value();
    s = s.play(2 * M + 1).value();
    return (s.z = 4);
}

// Test state.

// Regression test for https://www.codecup.nl/showgame.php?ga=164794
bool testInvalidMove() {
    Board b = 1152894691084558198;
    auto s = State{b, 30};
    s.play(20).value();
    s.play(21).value();
    s.play(32).value();
    s.play(37).value();
    s.play(43).value();

    auto moves = s.validMoves();
    return moves.size() == 5;
}

// Test mcts.

bool testSpeedSimulate() {
    auto s = State{};
    auto start = high_resolution_clock::now();
    ulong wins = 0;
    for (int i = 0; i < 1000; i++) {
        wins += s.simulate();
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
//    log(to_string(duration.count()));
    // 50ms, 50 microseconds per simulation
    return duration.count() < 50000 && wins > 100 && wins < 900;
}

bool testScore() {
    Node p = Node();
    p.n = 31;
    p.nLog = std::log(p.n);
    auto child1 = Node(&p, Move{});
    child1.n = 11;
    child1.nLog = std::log(child1.n);
    auto child2 = Node(&p, Move{});
    child2.n = 10;
    child2.nLog = std::log(child2.n);
    auto child3 = Node(&p, Move{});
    child3.n = 10;
    child3.nLog = std::log(child3.n);
    child3.w = 1; // 9 wins from parent
    return child1.score(2) < child2.score(2)
           && child2.score(2) > child3.score(2)
           && child1.score(1) < child1.score(2);
}


bool testBestMoveWithThreeSquares() {
    auto s = State{1152921504602650622, 0};
    Mcts mcts = Mcts(s, 2.0);
    auto node = mcts.mcts(10.0, 100);

    return node->m == 0 && node->winChances() > 0.9;
}

bool testBestMoveWithTwoSquares() {
    State s{1152921504606844926, 0};
    Mcts mcts = Mcts(s, 2.0);

    auto node = mcts.mcts(10.0, 100);

    // both moves are a win
    return node->winChances() > 0.9;
}

bool testLosingWithThreeSquares() {
    State s{1152921504602650622, 1 << 3};
    Mcts mcts = Mcts(s, 2.0);

    auto node = mcts.mcts(10.0, 100);

    return node->winChances() < 0.1;
}

bool testWinAgainstRandomized() {
    State s{0, 0};
    bool is_mcts = true;
    Mcts mcts = Mcts(s, 2.0);

    optional<pair<Move, State>> moves;
    do {
        Move m;
        if (is_mcts) {
            m = mcts.mcts(100000.0, 10000)->m;
            s = mcts.updateRoot(m);
        } else {
            m = s.randomMove().value().first;
            s = mcts.updateRoot(m);
        }
        moves = s.randomMove();
        is_mcts = !is_mcts;
    } while (!moves.has_value());

    return !is_mcts;
}

[[maybe_unused]] void test() {
    std::cout << "Tests" << std::endl;
    std::cout << testPlayWall()
              << " testPlayWall" << std::endl;
    std::cout << testCannotPutWallTwiceOnSamePosition()
              << " testCannotPutWallTwiceOnSamePosition" << std::endl;
    std::cout << testPlaceAllBorderWalls()
              << " testPlaceAllBorderWalls" << std::endl;
    std::cout << testPlayingOneTile()
              << " testPlayingOneTile" << std::endl;
    std::cout << testPlayingSquareTile()
              << " testPlayingSquareTile" << std::endl;
    std::cout << testInvalidMove()
              << " testInvalidMove" << std::endl;
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

// 2021-01-17 8.35193s, 8.35193us per iter
// 2021-01-18 6.66694s, 6.66694us per iter - vector instead of unordered_map
// 2021-01-18 6.32424s, 6.32424us per iter - remember excluded zones
// 2021-01-19 6.12408s, 6.12408us per iter - use bitset._Find_next()
// 2021-01-19 5.97594s, 5.97594us per iter - don't save states in tree
// 2021-01-21 4.46186s, 4.46186us per iter - only add 1 node during expand
[[maybe_unused]] void profile() {
    const uint iter = 1000000;
    State s{0, 0};
    Mcts mcts = Mcts(s, 2.0);

    auto start = high_resolution_clock::now();
    mcts.mcts(100.0, iter);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start).count();
    std::cerr << duration / 1000.0 / 1000.0 << "s, "
              << ((float) duration / (float) iter) << "us per iter"
              << std::endl;
}

int main() {
//    test();
    game(0.5, 1);
//    profile();

    return 0;
}
