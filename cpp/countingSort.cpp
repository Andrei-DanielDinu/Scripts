#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>  // For table formatting

using namespace std;

void countingSort(vector<int>& arr) {

}

int main() {
    cout<<"Welcome to SelectionSort algorithm"<<endl<<endl;

    vector<int> arr = {4, 2, 2, 8, 3, 3, 1};
    int n = arr.size();
    cout<<"Total count of elements is "<<n<<endl<<endl;

    cout << "Unsorted array: "<<endl;
    for(int i=0;i<n;i++)
        cout<<"On position "<<i<<" is the value "<<arr[i]<<"; arr["<<i<<"]="<<arr[i]<<endl;

    cout<<endl;

    cout<<"In the Counting Sort algorithm we:"<<endl<<"1. Store the frequency of each element of the array;"<<endl<<"2. Use cumulative counts for a simple indication of the numbers that are repeated or missing from the array."<<endl<<"3. Can determine the position and the frequency of the element in ascending order based on the cumulative count."<<endl<<endl;

    // Table Header
    cout << left 
              << setw(20) << "Complexity"
              << setw(15) << "Best Case"
              << setw(15) << "Average Case"
              << setw(15) << "Worst Case"
              << "\n";

    // Separator Line
    cout << setfill('-') << setw(65) << "" << setfill(' ') << "\n";

    // Time Complexity Row
    cout << left
              << setw(20) << "Time Complexity"
              << setw(15) << "O(n+k)"
              << setw(15) << "O(n+k)"
              << setw(15) << "O(n+k)"
              << "\n";

    // Space Complexity Row
    cout << setw(20) << "Space Complexity"
              << setw(15) << "O(n+k)"
              << setw(15) << "O(n+k)"
              << setw(15) << "O(n+k)"
              << "\n";
    
    cout<<endl;
    
    // Find max and min values
    int max_val = *max_element(arr.begin(), arr.end());
    cout<<"Maximum value found in array is "<<max_val<<endl;
    int min_val = *min_element(arr.begin(), arr.end());
    cout<<"Minimum value found in array is "<<min_val<<endl;
    int range = max_val - min_val + 1;
    cout<<"Array range is "<<range<<endl;

    // Initialize count and output arrays
    vector<int> count(range), output(arr.size());

    cout<<"Frequency array of the elements:"<<endl;
    // Count frequencies
    for (int num : arr) 
        count[num - min_val]++;
    for (int i = 0; i < range; i++)
        cout<<count[i]<<" ";
    cout<<endl<<endl;

    cout<<"Cumulative count array of the elements:"<<endl;
    // Compute prefix sums (cumulative counts)
    for (int i = 1; i < range; i++)
        count[i] += count[i - 1];
    for (int i = 0; i < range; i++)
        cout<<count[i]<<" ";
    cout<<endl<<endl;

    cout << "Creating output array: "<<endl;
    // Build the output array
    for (int i = arr.size() - 1; i >= 0; i--) {
        cout<<"We take the element from position "<<i<<" with value "<<arr[i]<<endl;
        cout<<"We calculate the distance in reference with the element with minimum value "<<arr[i]<<" - "<<min_val<<" = "<<arr[i] - min_val<<endl;
        cout<<"The appropriate position of the element "<<arr[i]<<" is "<<count[arr[i] - min_val] - 1<<endl;
        output[count[arr[i] - min_val] - 1] = arr[i];
        cout<<"We decrease the correspondent cumulative count "<<count[arr[i] - min_val]<<" - "<<"1"<<endl;
        count[arr[i] - min_val]--;
        cout<<endl;
        for(int i=0;i<n;i++)
            cout<<output[i]<<" ";

        cout<<endl;
        cout<<endl;
    }

    return 0;
}