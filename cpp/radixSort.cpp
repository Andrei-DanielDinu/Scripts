#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

// Counting Sort subroutine for Radix Sort
void countingSort(vector<int>& arr, int exp, bool debug = true) {
    vector<int> output(arr.size());
    vector<int> count(10, 0);  // Digits 0-9

    if (debug) {
        cout << "\n Processing digit place: " << exp
             << " (1s, 10s, 100s, etc.)\n";
    }

    // 1. Count frequencies of current digit
    for (int num : arr) {
        int digit = (num / exp) % 10;
        count[digit]++;
        if (debug) {
            cout << "  Number " << setw(3) << num << " -> Digit " << digit
                 << " | Count[" << digit << "] = " << count[digit] << endl;
        }
    }

    if (debug) {
        cout << "\n  Frequency array: [ ";
        for (int i = 0; i < 10; i++) cout << count[i] << " ";
        cout << "]\n";
    }

    // 2. Compute cumulative counts (prefix sum)
    for (int i = 1; i < 10; i++) {
        count[i] += count[i - 1];
    }

    if (debug) {
        cout << "  Cumulative counts: [ ";
        for (int i = 0; i < 10; i++) cout << count[i] << " ";
        cout << "]\n";
    }

    // 3. Build output array (stable sort)
    if (debug) cout << "\n  Placing numbers in sorted order:\n";
    for (int i = arr.size() - 1; i >= 0; i--) {
        int digit = (arr[i] / exp) % 10;
        output[count[digit] - 1] = arr[i];
        if (debug) {
            cout << "    Moved " << setw(3) << arr[i] << " to position "
                 << (count[digit] - 1) << endl;
        }
        count[digit]--;
    }

    arr = output;

    if (debug) {
        cout << "\n  After sorting this digit: [ ";
        for (int num : arr) cout << num << " ";
        cout << "]\n";
        cout << "-----------------------------------------------\n";
    }
}

void radixSort(vector<int>& arr, bool debug = true) {
    if (arr.empty()) return;

    if (debug) {
        cout << "Initial array: [ ";
        for (int num : arr) cout << num << " ";
        cout << "]\n\n";
    }

    int max_val = *max_element(arr.begin(), arr.end());
    if (debug) {
        cout << "Max value: " << max_val << " (will process "
             << to_string(max_val).length() << " digits)\n";
    }

    // Process each digit from LSD to MSD
    for (int exp = 1; max_val / exp > 0; exp *= 10) {
        countingSort(arr, exp, debug);
    }
}

int main() {
    vector<int> arr = {170, 45, 75, 90, 802, 24, 2, 66};

    // Complexity Table
    cout << "|-----------------------|--------------|--------------|-----------"
            "---|\n";
    cout << "│ " << left << setw(22) << "Complexity"
         << "│ " << setw(13) << "Best Case"
         << "│ " << setw(13) << "Average Case"
         << "│ " << setw(13) << "Worst Case" << "│\n";
    cout << "|-----------------------|--------------|--------------|-----------"
            "---|\n";
    cout << "│ " << setw(22) << "Time (Radix Sort)"
         << "│ " << setw(13) << "O(d(n+b))"
         << "│ " << setw(13) << "O(d(n+b))"
         << "│ " << setw(13) << "O(d(n+b))" << "│\n";
    cout << "│ " << setw(22) << "Space"
         << "│ " << setw(13) << "O(n+b)"
         << "│ " << setw(13) << "O(n+b)"
         << "│ " << setw(13) << "O(n+b)" << "│\n";
    cout << "|-----------------------|--------------|--------------|-----------"
            "---|\n\n";

    // Sort with debug output
    radixSort(arr, true);

    cout << "\nFinal sorted array: [ ";
    for (int num : arr) cout << num << " ";
    cout << "]\n";

    return 0;
}