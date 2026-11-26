// Independent small-weight control. Enumerate all matchings between distinct
// alternating cycles, then test the other two colours only on completed graphs.
// Usage: check_catalogues WEIGHT OUTPUT. Practical through weight 22.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

class Control {
    using Key = std::vector<unsigned char>;
    int n;
    std::vector<std::array<int,3>> graph;
    std::vector<int> cycle;
    std::set<Key> seen;

    bool valid(int omitted) const {
        std::vector<int> label(n, -1);
        for (int start = 0; start < n; ++start) {
            if (label[start] >= 0) continue;
            int v = start, c = (omitted+1)%3;
            do {
                label[v] = start;
                v = graph[v][c];
                c = 3-omitted-c;
            } while (v != start);
        }
        for (int v = 0; v < n; ++v)
            if (label[v] == label[graph[v][omitted]]) return false;
        return true;
    }
    void record() {
        if (!valid(0) || !valid(1)) return;
        Key minimum;
        for (int start = 0; start<n; ++start) {
            std::array<int,32> order{}, id{};
            id.fill(-1); id[start] = 0; order[0] = start;
            int filled = 1;
            Key key;
            for (int i = 0; i<filled; ++i) for (int v : graph[order[i]]) {
                if (id[v] < 0) { id[v] = filled; order[filled++] = v; }
                key.push_back(static_cast<unsigned char>(id[v]));
            }
            if (filled != n) return;
            if (minimum.empty() || key < minimum) minimum = std::move(key);
        }
        seen.insert(std::move(minimum));
    }
    void match(std::uint32_t open) {
        if (!open) { record(); return; }
        int a = __builtin_ctz(open);
        open ^= 1u << a;
        for (auto rest = open; rest; rest &= rest-1) {
            int b = __builtin_ctz(rest);
            if (cycle[a] == cycle[b]) continue;
            graph[a][2] = b; graph[b][2] = a;
            match(open ^ (1u << b));
        }
    }
    void partitions(int remaining, int least, std::vector<int>& sizes) {
        if (remaining) {
            for (int len=least; len<=remaining && len<=n/2; len+=2) {
                sizes.push_back(len);
                partitions(remaining-len, len, sizes);
                sizes.pop_back();
            }
            return;
        }
        int pos=0;
        for (int len:sizes) {
            for (int j=0; j<len; ++j) {
                int a=pos+j, b=pos+(j+1)%len;
                graph[a][j%2]=b; graph[b][j%2]=a; cycle[a]=pos;
            }
            pos+=len;
        }
        match((1u<<n)-1);
    }
public:
    explicit Control(int w):n(w),graph(w),cycle(w) {}
    void output(std::ostream& out) {
        std::vector<int> sizes;
        partitions(n,4,sizes);
        out << seen.size() << '\n';
        for (const auto& key:seen) {
            out << n;
            for (int x:key) out << ' ' << x;
            out << '\n';
        }
        if (!out) throw std::runtime_error("Could not write control catalogue");
    }
};
int main(int argc,char** argv) {
    if (argc != 3) return 2;
    int n = std::stoi(argv[1]);
    if (n < 2 || n > 22 || n%2) return 2;
    std::ofstream out(argv[2]);
    Control(n).output(out);
}
