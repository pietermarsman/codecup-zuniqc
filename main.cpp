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
class State {
public:
    Board b = 0;
    Zones z = 0;

    // State functions.
    [[nodiscard]] optional<State> play(Move m) const {
        const auto maybe_b = put(b, m);
        if (!maybe_b.has_value()) {
            return nullopt;
        }

        Board bb = fillEnclosed(maybe_b.value());
        uint zone = countAddedZones(b, bb);
        if ((z >> zone & 1) == 1) {
            return nullopt;
        }
        Zones zz = updateZones(z, zone);
        return State{bb, zz};
    }

    [[nodiscard]] unordered_map<Move, State> validMoves() const {
        unordered_map<Move, State> moves;
        bitset<N> bits = b;
        for (int i = 0; i < N; ++i) {
            if (bits[i] == 0) {
                auto maybe_s = play(i);
                if (maybe_s.has_value()) {
                    moves[i] = *maybe_s;
                }
            }
        }
        return moves;
    }

    [[nodiscard]] optional<pair<Move, State>> randomMove() const {
        bitset<N> bBits = b;
        uint shift = random() % N;
        for (uint i = 0; i < N; i++) {
            uint j = (i + shift) % N;
            if (bBits[j] == 0) {
                auto maybe_s = play(j);
                if (maybe_s.has_value()) {
                    const auto ms = pair<Move, State>{j, maybe_s.value()};
                    return ms;
                }
            }
        }
        return nullopt;
    }

    [[nodiscard]] list<uint> listZones() const {
        list<uint> zoneSizes;
        auto zBit = (bitset<64>) z;
        for (ulong i = 0; i < zBit.size(); i++) {
            if (zBit[i]) {
                zoneSizes.push_back(i);
            }
        }
        return zoneSizes;
    }

    ~State() = default;

    // Required for hash map.
    bool operator==(const State &s) const {
        return b == s.b && z == s.z;
    }

    [[nodiscard]] ulong simulate() const {
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

    static Board proceed(Board b) {
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

    static Board fillEnclosed(Board b) {
        Board bb = (b ^ FULL) & BORDER;
        Board previous;
        do {
            previous = bb;
            bb = proceed(bb) & (b ^ FULL);
        } while (bb != previous);

        return b | (bb ^ FULL);
    }

    static uint countOnes(Board b) {
        return ((bitset<64>) b).count();
    }

    // Zone functions.
    static uint countEnclosedTiles(Board b) {
        b = fillEnclosed(b);
        Board shiftBack = b & (b >> M) & (b >> W) & (b >> (W + 1));
        return countOnes(b & shiftBack & HORIZONTAL & ~BOTTOM);
    }

    static uint countAddedZones(Board b, Board bb) {
        return countEnclosedTiles(bb) - countEnclosedTiles(b);
    }

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
    State s = State{};
    uint depth = 0;

    bool expanded = false;
    vector<shared_ptr<Node>> children;

    ulong n = 0;
    ulong w = 0;

    // Create root node of t; no p
    explicit Node(State state) {
        s = state;
    }

    // Create intermediate node of t
    explicit Node(Node *parent, Move move, State state) {
        p = parent;
        m = move;
        s = state;
        depth = parent->depth + 1;
    }

    // mcts functions.
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

    // compute score from perspective of parent
    [[nodiscard]] float score(float c, float eps = 0.0001) const {
        float n_eps = n + eps;
        float n_parent_eps = p->n + eps;
        float wins = (float) n - w;

        return ((float) wins / n_eps)
               + c * sqrtf(std::log(n_parent_eps) / n_eps);
    }

    // wins from perspective of parent
    [[nodiscard]] ulong wins() const {
        return n - w;
    }

    void expand() {
        assert (!expanded);
        for (auto move_state : s.validMoves()) {
            auto child = make_shared<Node>(this, move_state.first,
                                           move_state.second);
            children.push_back(child);
        }
        expanded = true;
    }

    void backPropagate(ulong win) {
        n++;
        w += win;
        if (p) {
            p->backPropagate(1 - win);
        }
    }
};


class Mcts {
public:
    explicit Mcts(State s, float constant) {
        t = make_shared<Node>(s);
        c = constant;
    }

    // Execute mcts searches for fixed duration.
    shared_ptr<Node>
    mcts(float duration_secs, float max_iterations = INFINITY) {
        auto start = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(
                high_resolution_clock::now() - start);

        ulong i = 0;
        do {
            i++;
            search(t);
            duration = duration_cast<microseconds>(
                    high_resolution_clock::now() - start);
        } while ((duration.count() / 1000.0 / 1000.0) < duration_secs
                 && i < max_iterations);

        return getBestMove();
    }

    void updateRoot(Move move) {
        shared_ptr<Node> node;
        for (const auto &child: t->children) {
            if (child->m == move) {
                node = child;
            }
        }
        t = node;
        t->p = nullptr;
    }

    float timeBudget(float total_time) {
        float time_left = 29.0F - total_time;
        float moves_left = t->s.validMoves().size();
        float max_move_time = time_left / 2.0F;
        float min_move_time = time_left / (moves_left + 1);

        if (moves_left < 20) {
            return min(min_move_time * 4.0F, max_move_time);
        } else if (moves_left < 35) {
            return min(min_move_time * 2.0F, max_move_time);
        } else {
            return min(min_move_time * 0.1F, max_move_time);
        }
    }

private:
    shared_ptr<Node> t;
    float c;

    void search(shared_ptr<Node> n) const {
        while (n->expanded && !n->children.empty()) {
            n = n->bestMove(c);
        }
        if (!n->expanded) {
            // stopped because it is not yet expanded, not because it has no
            // children.
            n->expand();
        }
        ulong win = n->s.simulate();
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
    float win_chances = (float) node->wins() / (float) node->n;
    log("Chances of winning are " + to_string(node->wins()) + " / "
        + to_string(node->n) + " = " + to_string(win_chances) + "\n");
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
        list<uint> zoneSizes = s.listZones();
        log("Enclosed zones: ");
        for (uint i: zoneSizes) {
            log(to_string(i) + ", ");
        }
        log("\n");
    }
    if (verbose >= 0) {
        unordered_map<Move, State> moves = s.validMoves();
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

void game() {
    float total_ms = 0;
    Mcts mcts = Mcts(State{}, sqrt(2.0F));

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
        logMcts(best_node);
        logState(best_node->s, 1);
        mcts.updateRoot(best_node->m);
        write(INT_TO_POSITION[best_node->m]);

        auto end_time = high_resolution_clock::now();
        total_ms += duration_cast<milliseconds>(end_time - start_time).count();
        logState(best_node->s);
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
    State s1 = s.play(20).value();
    State s2 = s.play(21).value();
    State s3 = s.play(32).value();
    State s4 = s.play(37).value();
    State s5 = s.play(43).value();

    unordered_map<Move, State> moves = s.validMoves();
    return moves.size() == 5
           && moves[20] == s1
           && moves[21] == s2
           && moves[32] == s3
           && moves[37] == s4
           && moves[43] == s5;
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
    Node p = Node(State{});
    p.n = 31;
    auto child1 = Node(&p, Move{}, State{});
    child1.n = 11;
    auto child2 = Node(&p, Move{}, State{});
    child2.n = 10;
    auto child3 = Node(&p, Move{}, State{});
    child3.n = 10;
    child3.w = 1; // 9 wins from parent
    return child1.score(2) < child2.score(2)
           && child2.score(2) > child3.score(2)
           && child1.score(1) < child1.score(2);
}


bool testBestMoveWithThreeSquares() {
    auto s = State{1152921504602650622, 0};
    Mcts mcts = Mcts(s, 2.0);
    auto node = mcts.mcts(10.0, 100);
    float chances = (float) node->wins() / (float) node->n;

    return node->m == 0 && chances > 0.9;
}

bool testBestMoveWithTwoSquares() {
    State s{1152921504606844926, 0};
    Mcts mcts = Mcts(s, 2.0);

    auto node = mcts.mcts(10.0, 100);
    float chances = (float) node->wins() / (float) node->n;

    // both moves are a win
    return chances > 0.9;
}

bool testLosingWithThreeSquares() {
    State s{1152921504602650622, 1 << 3};
    Mcts mcts = Mcts(s, 2.0);

    auto node = mcts.mcts(10.0, 100);
    float chances = (float) node->wins() / (float) node->n;

    return chances < 0.1;
}

bool testWinAgainstRandomized() {
    State s{0, 0};
    bool is_mcts = true;
    Mcts mcts = Mcts(s, 2.0);

    unordered_map<Move, State> moves;
    do {
        Move m;
        if (is_mcts) {
            m = mcts.mcts(10.0, 100)->m;
            mcts.updateRoot(m);
        } else {
            m = s.randomMove().value().first;
            mcts.updateRoot(m);
        }
        s = s.play(m).value();
        moves = s.validMoves();
        is_mcts = !is_mcts;
    } while (!moves.empty());

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

int main() {
    test();
//    game();

    // todo figure out good simulate() for start of game

    return 0;
}
