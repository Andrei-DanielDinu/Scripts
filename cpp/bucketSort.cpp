#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

void bucketSort(vector<float>& arr, bool debug = true) {
    if (arr.empty()) return;

    // 1. Initialize buckets
    int n = arr.size();
    vector<vector<float>> buckets(n);

    if (debug) {
        cout << "Initial array: [ ";
        for (float num : arr) cout << num << " ";
        cout << "]\n\n";
        cout << "Step 1: Scatter elements into buckets\n";
    }

    // 2. Scatter elements into buckets
    for (float num : arr) {
        int bucketIdx = n * num;  // For floats in [0, 1)
        buckets[bucketIdx].push_back(num);

        if (debug) {
            cout << "  Value " << fixed << setprecision(2) << num
                 << " -> Bucket " << bucketIdx << "\n";
        }
    }

    if (debug) {
        cout << "\nStep 2: Sort individual buckets\n";
        for (int i = 0; i < n; i++) {
            if (!buckets[i].empty()) {
                cout << "  Bucket " << i << ": [ ";
                for (float num : buckets[i]) cout << num << " ";
                cout << "] -> Sorted: [ ";
                sort(buckets[i].begin(), buckets[i].end());
                for (float num : buckets[i]) cout << num << " ";
                cout << "]\n";
            }
        }
    } else {
        for (auto& bucket : buckets) {
            sort(bucket.begin(), bucket.end());
        }
    }

    // 3. Gather elements from buckets
    int idx = 0;
    if (debug) cout << "\nStep 3: Concatenate buckets\n";
    for (const auto& bucket : buckets) {
        for (float num : bucket) {
            arr[idx++] = num;
            if (debug)
                cout << "  Placing " << num << " -> Position " << idx - 1
                     << "\n";
        }
    }
}

int main() {
    vector<float> arr = {0.42, 0.32, 0.75, 0.12, 0.89, 0.63};

    // Complexity Table
    cout << "|-----------------------|--------------|--------------|-----------"
            "---|\n";
    cout << "| " << left << setw(22) << "Complexity"
         << "| " << setw(13) << "Best Case"
         << "| " << setw(13) << "Average Case"
         << "| " << setw(13) << "Worst Case" << "|\n";
    cout << "|-----------------------|--------------|--------------|-----------"
            "---|\n";
    cout << "| " << setw(22) << "Time (Bucket Sort)"
         << "| " << setw(13) << "O(n + k)"
         << "| " << setw(13) << "O(n + k)"
         << "| " << setw(13) << "O(n*n)" << "|\n";
    cout << "| " << setw(22) << "Space"
         << "| " << setw(13) << "O(n + k)"
         << "| " << setw(13) << "O(n + k)"
         << "| " << setw(13) << "O(n + k)" << "|\n";
    cout << "|-----------------------|--------------|--------------|-----------"
            "---|\n\n";

    // Sort with debug output
    bucketSort(arr, true);

    cout << "\nFinal sorted array: [ ";
    for (float num : arr) cout << num << " ";
    cout << "]\n";

    return 0;
}