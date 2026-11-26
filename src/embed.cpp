// Exhaustively embed each support pattern in one specified CSS sector.
// Input from prepare.py: checks of each qubit, base coordinates, stabilizer rows.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

using Triple = std::array<int, 3>;
using Bits = std::vector<std::uint64_t>;
constexpr int lift = 32, side = 8;
struct Step {
    int vertex, parent, colour, coordinate_source;
    std::array<int, 2> close_at{}, close_colour{};
    int closures = 0;
};

class Search {
    int n, m;
    std::vector<Triple> check, coordinate;
    std::vector<std::array<int, side>> line;
    std::vector<Bits> basis;
    std::vector<bool> used;
    struct Node {
        Step step{};
        int child = -1, sibling = -1, pattern = -1;
    };
    std::vector<Node> tree;
    std::vector<Triple> stabilizers;
    std::vector<unsigned char> occupancy;
    int target_weight = 0, found_pattern = -1;
    bool minimal_only = false;
    static bool equivalent(const Step& a, const Step& b) {
        return a.parent == b.parent && a.colour == b.colour &&
               a.coordinate_source == b.coordinate_source && a.closures == b.closures &&
               a.close_at == b.close_at && a.close_colour == b.close_colour;
    }
    void add(const std::vector<Step>& order, int pattern) {
        int node = 0;
        for (size_t i = 1; i < order.size(); ++i) {
            int child = tree[node].child;
            while (child >= 0 && !equivalent(tree[child].step, order[i])) child = tree[child].sibling;
            if (child < 0) {
                child = static_cast<int>(tree.size());
                int previous = tree[node].child;
                tree[node].child = child;
                tree.push_back({order[i], -1, previous, -1});
            }
            node = child;
        }
        if (tree[node].pattern != -1) throw std::runtime_error("Duplicate search plan");
        tree[node].pattern = pattern;
    }
    std::vector<int> assigned;
    std::uint64_t leaves = 0;
    std::vector<int> witness;

    int pivot(const Bits& v, int from) const {
        for (int i = from; i >= 0; --i)
            if (v[i]) return 64 * i + 63 - __builtin_clzll(v[i]);
        return -1;
    }
    bool reduce(Bits& word, bool insert) {
        int start = static_cast<int>(word.size()) - 1;
        for (;;) {
            int p = pivot(word, start);
            if (p < 0) return true;
            if (basis[p].empty()) {
                if (insert) basis[p] = std::move(word);
                return false;
            }
            start = p / 64;
            for (int i = 0; i <= start; ++i) word[i] ^= basis[p][i];
        }
    }
    static void toggle(Bits& bits, int q) {
        bits.at(q / 64) ^= std::uint64_t{1} << (q % 64);
    }
    // Choose only the visiting order here. Every eligible qubit is tried in visit().
    static std::vector<Step> plan(const std::vector<Triple>& adj) {
        int w = static_cast<int>(adj.size());
        std::vector<Triple> component(w, {-1,-1,-1});
        for (int c = 0; c < 3; ++c) for (int i = 0; i < w; ++i) {
            if (component[i][c] >= 0) continue;
            int v = i, colour = (c + 1) % 3;
            do {
                component[v][c] = i;
                v = adj[v][colour];
                colour = 3 - c - colour;
            } while (v != i);
        }
        double best = std::numeric_limits<double>::infinity();
        std::vector<Step> result;
        for (int root = 0; root < std::min(3,w); ++root) {
            std::vector<int> visited(w, -1);
            std::vector<Triple> first(w, {-1,-1,-1});
            std::vector<Step> order{{root,-1,-1,-1}};
            visited[root] = 0;
            for (int c = 0; c < 3; ++c) first[component[root][c]][c] = 0;
            double population = 1, cost = 0;
            while (static_cast<int>(order.size()) < w) {
                int best_priority = -100;
                Step chosen{};
                for (int v = 0; v < w; ++v) {
                    if (visited[v] >= 0) continue;
                    Step s{v,-1,-1,-1};
                    for (int c = 0; c < 3; ++c) {
                        int other = visited[adj[v][c]];
                        if (other < 0) continue;
                        if (s.parent < 0) {
                            s.parent = other; s.colour = c;
                            s.coordinate_source = first[component[v][c]][c];
                        } else {
                            s.close_at[s.closures] = other;
                            s.close_colour[s.closures++] = c;
                        }
                    }
                    if (s.parent < 0) continue;
                    int priority = 10*s.closures + 4*(s.coordinate_source >= 0);
                    if (priority > best_priority) { best_priority = priority; chosen = s; }
                }
                if (best_priority < 0) throw std::runtime_error("Disconnected pattern");
                population *= chosen.coordinate_source < 0 ? 7.0 : 1.0;
                cost += population;
                population *= (chosen.closures == 2 ? 0.0144 : chosen.closures == 1 ? 0.12 : 1.0);
                int at = static_cast<int>(order.size());
                visited[chosen.vertex] = at;
                for (int c = 0; c < 3; ++c) {
                    int& source = first[component[chosen.vertex][c]][c];
                    if (source < 0) source = at;
                }
                order.push_back(chosen);
            }
            if (cost < best) { best = cost; result = std::move(order); }
        }
        return result;
    }

    bool visit(int node, int depth) {
        if (depth == target_weight) {
            ++leaves;
            Bits word((n + 63)/64);
            for (int q : assigned) toggle(word, q);
            if (reduce(word, false)) return false;
            found_pattern = tree[node].pattern;
            witness = assigned;
            std::sort(witness.begin(), witness.end());
            return true;
        }
        for (int child = tree[node].child; child >= 0; child = tree[child].sibling) {
          const Step& s = tree[child].step;
        int r = check[assigned[s.parent]][s.colour];
        int begin = 0, end = side;
        if (s.coordinate_source >= 0) {
            begin = coordinate[assigned[s.coordinate_source] / lift][s.colour];
            end = begin + 1;
        }
        for (int t = begin; t < end; ++t) {
            int q = line[r][t];
            if (used[q]) continue;
            bool consistent = true;
            for (int i = 0; i < s.closures; ++i) {
                int c = s.close_colour[i];
                if (check[q][c] != check[assigned[s.close_at[i]]][c]) consistent = false;
            }
            if (!consistent) continue;
            if (minimal_only) {
                const auto& rows = stabilizers[q];
                if (occupancy[rows[0]] >= 4 || occupancy[rows[1]] >= 4 || occupancy[rows[2]] >= 4)
                    continue;
                for (int r : rows) ++occupancy[r];
            }
            assigned[depth] = q;
            used[q] = true;
            bool found = visit(child, depth + 1);
            if (minimal_only) for (int r : stabilizers[q]) --occupancy[r];
            used[q] = false;
            if (found) return true;
          }
        }
        return false;
    }
public:
    explicit Search(std::istream& in) {
        if (!(in >> n >> m) || n != 16384 || m != 6144) throw std::runtime_error("Bad dimensions");
        check.resize(n); coordinate.resize(n/lift); line.resize(m); basis.resize(n); used.resize(n);
        for (auto& row : line) row.fill(-1);
        for (auto& xyz : coordinate) for (int& value : xyz)
            if (!(in >> value) || value < 0 || value >= side) throw std::runtime_error("Invalid coordinate");
        for (int q = 0; q < n; ++q) for (int c = 0; c < 3; ++c) {
            int r;
            if (!(in >> r) || r < 0 || r >= m || r/(m/3) != c) throw std::runtime_error("Bad row or colour");
            check[q][c] = r;
            int& slot = line[r][coordinate[q/lift][c]];
            if (slot >= 0) throw std::runtime_error("Repeated coordinate in one check");
            slot = q;
        }
        for (const auto& row : line)
            if (std::count(row.begin(),row.end(),-1)) throw std::runtime_error("Incomplete base line");
        stabilizers.assign(n, {-1,-1,-1});
        occupancy.assign(m,0);
        for (int r = 0; r < m; ++r) {
            Bits word((n+63)/64);
            for (int i = 0; i < side; ++i) {
                int q;
                if (!(in >> q) || q < 0 || q >= n) throw std::runtime_error("Bad stabilizer");
                stabilizers[q][r/(m/3)] = r;
                toggle(word, q);
            }
            reduce(word, true);
        }
        int rank = std::count_if(basis.begin(),basis.end(),[](const Bits& b){ return !b.empty(); });
        if (rank != 6121) throw std::runtime_error("Stabilizer rank mismatch");
    }
    void search(std::istream& input, int weight, bool prune) {
        target_weight = weight; minimal_only = prune;
        int count;
        if (!(input >> count) || count < 0) throw std::runtime_error("Bad pattern count");
        constexpr int batch_size = 2048;
        for (int batch = 0; batch < count; batch += batch_size) {
            tree.clear();
            tree.emplace_back();
            int end = std::min(count, batch + batch_size);
            for (int index = batch; index < end; ++index) {
                int w;
                if (!(input >> w) || w != weight) throw std::runtime_error("Bad pattern weight");
                std::vector<Triple> adj(w);
                for (auto& links : adj) for (int& v : links)
                    if (!(input >> v) || v<0 || v>=w) throw std::runtime_error("Bad pattern edge");
                add(plan(adj), index);
            }
            assigned.resize(weight);
            for (int root = 0; root < n; root += lift) {
                assigned[0] = root; used[root] = true;
                if (prune) for (int r : stabilizers[root]) ++occupancy[r];
                bool found = visit(0, 1);
                if (prune) for (int r : stabilizers[root]) --occupancy[r];
                used[root] = false;
                if (found) {
                    std::cout << "{\"found\":true,\"completed_patterns\":" << batch
                              << ",\"witness_pattern_zero_based\":" << found_pattern
                              << ",\"embeddings_tested\":" << leaves << ",\"support\":[";
                    for (int i = 0; i < weight; ++i) std::cout << (i ? "," : "") << witness[i];
                    std::cout << "]}\n";
                    return;
                }
            }
        }
        std::string extra;
        if (input >> extra) throw std::runtime_error("Extra pattern rows");
        std::cout << "{\"found\":false,\"completed_patterns\":" << count
                  << ",\"embeddings_tested\":" << leaves << "}\n";
    }
};

int main(int argc, char** argv) {
    if (argc != 4 && argc != 5) return 2;
    if (argc == 5 && std::string(argv[4]) != "--minimal") return 2;
    std::ifstream matrix(argv[1]), patterns(argv[2]);
    Search solver(matrix);
    solver.search(patterns,std::stoi(argv[3]),argc == 5);
}
