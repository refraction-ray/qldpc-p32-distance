// Enumerate connected, coordinate-compatible, three-edge-coloured support patterns.
// Usage: patterns EVEN_WEIGHT OUTPUT. Vertices in the output are in coloured BFS order.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using Mask = std::uint32_t;
using Graph = std::vector<std::array<int, 3>>;
static Mask bit(int v) { return Mask{1} << v; }

// For a deleted colour c, vertices in any component using the other two
// colours must have the same base coordinate c. The endpoints of a c-edge
// must instead differ in that coordinate. Test this as components grow.
struct Components {
    std::array<int, 32> parent{};
    std::array<Mask, 32> members{}, forbidden{};
    struct Change { int a = -1, b = -1; Mask old_members = 0, old_forbidden = 0; };

    void init(const Graph& g, int deleted) {
        std::iota(parent.begin(), parent.end(), 0);
        for (int v = 0; v < static_cast<int>(g.size()); ++v) {
            int u = g[v][1 - deleted];
            int root = std::min(v, u);
            parent[v] = root;
            members[root] = bit(v) | bit(u);
            forbidden[root] = bit(g[v][deleted]) | bit(g[u][deleted]);
        }
    }
    int root(int v) const {
        while (parent[v] != v) v = parent[v];
        return v;
    }
    bool join(int x, int y, Change& change) {
        int a = root(x), b = root(y);
        if (a == b) return true;
        if (members[a] & forbidden[b]) return false;
        change = {a, b, members[a], forbidden[a]};
        parent[b] = a;
        members[a] |= members[b];
        forbidden[a] |= forbidden[b];
        return true;
    }
    void undo(Change change) {
        if (change.a < 0) return;
        parent[change.b] = change.b;
        members[change.a] = change.old_members;
        forbidden[change.a] = change.old_forbidden;
    }
};

class Catalog {
    using Key = std::vector<unsigned char>;
    int weight;
    Mask full;
    Graph g;
    std::vector<Mask> cycle;
    std::array<Components, 2> components;
    std::set<Key> keys;

    void record() {
        Key best;
        for (int root = 0; root < weight; ++root) {
            std::array<int, 32> id{}, order{};
            id.fill(-1);
            id[root] = 0;
            order[0] = root;
            int size = 1;
            Key key;
            for (int i = 0; i < size; ++i) {
                for (int next : g[order[i]]) {
                    if (id[next] < 0) { id[next] = size; order[size++] = next; }
                    key.push_back(static_cast<unsigned char>(id[next]));
                }
            }
            if (size != weight) return;
            if (best.empty() || key < best) best = std::move(key);
        }
        keys.insert(std::move(best));
    }

    void match(Mask open) {
        if (!open) { record(); return; }
        int a = __builtin_ctz(open);
        Mask candidates = open & ~cycle[a];
        Mask untouched_sizes = 0;
        while (candidates) {
            int b = __builtin_ctz(candidates);
            candidates &= candidates - 1;
            // A pristine alternating cycle is vertex-transitive under colour-
            // preserving maps. Equal-length pristine cycles may also be exchanged.
            if ((cycle[b] & open) == cycle[b]) {
                Mask size = bit(__builtin_popcount(cycle[b]));
                if (untouched_sizes & size) continue;
                untouched_sizes |= size;
            }
            Components::Change first, second;
            if (!components[0].join(a, b, first)) continue;
            if (components[1].join(a, b, second)) {
                g[a][2] = b;
                g[b][2] = a;
                match(open ^ bit(a) ^ bit(b));
                components[1].undo(second);
            }
            components[0].undo(first);
        }
    }

    void partition(int start, int minimum) {
        if (start == weight) {
            for (int c = 0; c < 2; ++c) components[c].init(g, c);
            match(full);
            return;
        }
        // Every colour-2 partner lies outside its own alternating 0/1 cycle.
        for (int len = minimum; len <= std::min(weight - start, weight / 2); len += 2) {
            Mask vertices = (bit(len) - 1) << start;
            for (int i = 0; i < len; ++i) {
                int a = start + i, b = start + (i + 1) % len;
                g[a][i % 2] = b;
                g[b][i % 2] = a;
                cycle[a] = vertices;
            }
            partition(start + len, len);
        }
    }
public:
    explicit Catalog(int n) : weight(n), full(bit(n) - 1), g(n), cycle(n) {}
    void write(std::ostream& out) {
        partition(0, 4);
        out << keys.size() << '\n';
        for (const auto& key : keys) {
            out << weight;
            for (int v : key) out << ' ' << v;
            out << '\n';
        }
        if (!out) throw std::runtime_error("Could not write the complete pattern list");
    }
};

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    int weight = std::stoi(argv[1]);
    if (weight < 2 || weight > 30 || weight % 2) return 2;
    std::ofstream out(argv[2]);
    Catalog(weight).write(out);
}
