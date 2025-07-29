#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

// Function to visualize current search range
void visualizeSearch(const vector<int>& arr, int low, int high, int pos,
                     int target) {
    int width = 70;  // Width of the ASCII plot
    int n = arr.size();

    // Determine value range
    int min_val = *min_element(arr.begin(), arr.end());
    int max_val = *max_element(arr.begin(), arr.end());
    double scale = static_cast<double>(width) / (max_val - min_val);

    cout << "\nCurrent Search Range: [" << arr[low] << " - " << arr[high]
         << "]\n";
    cout << "Estimated Position: index " << pos << " (value " << arr[pos]
         << ")\n";

    // Create ASCII plot
    for (int i = low; i <= high; i++) {
        // Calculate position in plot
        int position = static_cast<int>((arr[i] - min_val) * scale);

        cout << setw(4) << arr[i] << " |";
        for (int j = 0; j < width; j++) {
            if (j == position) {
                if (i == pos) {
                    cout << "x";  // Estimated position
                } else if (arr[i] == target) {
                    cout << "*";  // Target
                } else {
                    cout << "o";  // Data point
                }
            } else {
                cout << " ";
            }
        }
        cout << "|\n";
    }

    // Print position markers
    cout << "     +";
    for (int i = 0; i < width; i++) cout << "-";
    cout << "+\n";
    cout << "      ";
    for (int i = 0; i < width; i += 10) {
        cout << setw(10) << min_val + static_cast<int>(i / scale);
    }
    cout << "\n";
}

// Interpolation Search with visualization
int interpolationSearch(vector<int>& arr, int target, bool debug = true) {
    int low = 0;
    int high = arr.size() - 1;
    int steps = 0;

    if (debug) {
        cout << "|-------------------------------------------------------------"
                "--|\n";
        cout << "|                  INTERPOLATION SEARCH VISUALIZER            "
                " |\n";
        cout << "|-------------------------------------------------------------"
                "--|\n";
        cout << "| Target: " << setw(4) << target
             << "                                       |\n";
        cout << "| Array: ";
        for (int num : arr) cout << num << " ";
        cout << "       |\n";
        cout << "|-------------------------------------------------------------"
                "--|\n";
    }

    while (low <= high && target >= arr[low] && target <= arr[high]) {
        steps++;
        if (arr[high] == arr[low]) {
            if (arr[low] == target) return low;
            return -1;
        }

        // Calculate position using interpolation formula
        int pos = low + (((double)(target - arr[low]) * (high - low)) /
                         (arr[high] - arr[low]));
        pos = max(low, min(pos, high));  // Clamp to range

        if (debug) {
            cout << "\n\n=================== STEP " << steps
                 << " ===================\n";
            cout << "Low: " << setw(2) << low << " (value: " << arr[low]
                 << ")\n";
            cout << "High: " << setw(2) << high << " (value: " << arr[high]
                 << ")\n";
            cout << "Formula: pos = " << low << " + ((" << target << " - "
                 << arr[low] << ") * (" << high << " - " << low << ")) / ("
                 << arr[high] << " - " << arr[low] << ") = " << pos << "\n";

            visualizeSearch(arr, low, high, pos, target);
        }

        if (arr[pos] == target) {
            if (debug) {
                cout << "\nX TARGET FOUND AT INDEX " << pos << "! X\n";
                cout << "Steps taken: " << steps << "\n";
            }
            return pos;
        }

        if (arr[pos] < target) {
            low = pos + 1;
        } else {
            high = pos - 1;
        }
    }

    if (debug) cout << "\nX TARGET NOT FOUND! X\n";
    return -1;
}

int main() {
    // Uniformly distributed data (best case for Interpolation Search)
    vector<int> uniform_data = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    int target1 = 60;

    // Non-uniform data (worst case)
    vector<int> non_uniform_data = {1,   2,    100,  200,  300,
                                    500, 1000, 2000, 5000, 10000};
    int target2 = 500;

    // Complexity Table
    cout << "|---------------------------|---------------------------|\n";
    cout << "|      Time Complexity      |      Space Complexity     |\n";
    cout << "|-------------|-------------|-------------|-------------|\n";
    cout << "|   Best      |  Average    |   Worst     |   All Cases |\n";
    cout << "|-------------|-------------|-------------|-------------|\n";
    cout << "|   O(1)      | O(log log n)|   O(n)      |    O(1)     |\n";
    cout << "|-------------|-------------|-------------|-------------|\n\n";

    cout << "\n================ UNIFORM DISTRIBUTION (BEST CASE) "
            "================\n";
    int result1 = interpolationSearch(uniform_data, target1, true);

    cout << "\n\n============== NON-UNIFORM DISTRIBUTION (WORST CASE) "
            "==============\n";
    int result2 = interpolationSearch(non_uniform_data, target2, true);

    return 0;
}