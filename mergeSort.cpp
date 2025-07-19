#include <iostream>
#include <iomanip>  // For table formatting
#include <vector>

using namespace std;

void merge(vector<int>& arr, int left, int mid, int right) {
    //we put together what once was whole.
    int n1 = mid - left + 1;
    cout<<endl<<"First segment is "<<n1<<" elements long";
    int n2 = right - mid;
    cout<<endl<<"Second segment is "<<n1<<" elements long";

    vector<int> L(n1), R(n2);
    
    cout<<endl<<"First segment has elements: ";
    for (int i = 0; i < n1; i++)
        {L[i] = arr[left + i];
        cout<<L[i]<<" ";}
    cout<<endl<<L[n1-1]<<" is mid"<<endl;

    cout<<endl<<"Second segment has elements: ";
    for (int j = 0; j < n2; j++)
        {R[j] = arr[mid + 1 + j];
        cout<<R[j]<<" ";}

    int i = 0, j = 0, k = left;
    //we are so back(literally)
    cout<<endl;
    //now basically we cherrypick whatever we like from the other's basket.
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) { //Sometimes we're looking for something small
            arr[k] = L[i]; //where we expect it to be
            cout<<"Array["<<k<<"]="<<L[i]<<" found on left side."<<endl;
            i++;//rather disoriented than lost
        } else {
            arr[k] = R[j]; //but also where it shouldn't be 
            cout<<"Array["<<k<<"]="<<R[j]<<" found on right side."<<endl;
            j++; //keep searching for lost small numbers
        }
        cout<<arr[k]<<" is moved on Position "<<k<<" for the moment."<<endl;
        k++; //I found it, carry on
    
    }
    cout<<endl<<"Add what's left: + ";
    while (i < n1) {//complete with what's left
        arr[k] = L[i];
        cout<<arr[k]<<" on position "<<k<<" , ";
        i++;
        k++;
    }

    cout<<endl<<"Add what's right: + ";
    while (j < n2) {//and fill it right
        arr[k] = R[j];
        cout<<arr[k]<<" on position "<<k<<" , ";
        j++;
        k++;
    }
    cout<<endl<<endl;
}

void DivideEtImpera(vector<int>& arr, int left, int right) {
    if (left < right) {
        //If we still have more than one number to look at
        int mid = left + (right - left) / 2;
        cout<<"Calculate the mid between "<<left<<" and "<<right<<endl<<mid<<" = "<<left<<" + ("<<right<<" - "<<left<<") / 2"<<endl;
        //we find the middle
        DivideEtImpera(arr, left, mid);
        //mid is the new right
        DivideEtImpera(arr, mid + 1, right);
        //but it appears as left from the "true" right perspective
        merge(arr, left, mid, right);
        //after the halves sort their matters, we conciliate the parts (who gets the house and who gets the dog)
    }
    else
        cout<<endl<<"Ups, our division went too far: left Position = right Position = "<<left<<endl<<endl;
}

int main()
{
    cout<<"Welcome to MergeSort algorithm"<<endl<<endl;

    vector<int> arr = {64, 34, 25, 12, 22, 11, 90};
    int n = arr.size();
    cout<<"Total count of elements is "<<n<<endl<<endl;

    cout << "Unsorted array: "<<endl;
    for(int i=0;i<n;i++)
        cout<<"On position "<<i<<" is the value "<<arr[i]<<"; arr["<<i<<"]="<<arr[i]<<endl;

    cout<<endl;

    cout<<"In the Merge Sort algorithm we break down the problem into smaller pieces and then we compare adjacent halves."<<endl<<"After that, we complete with what's left."<<endl<<endl;

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
              << setw(15) << "O(n*log(n))"
              << setw(15) << "O(n*log(n))"
              << setw(15) << "O(n*log(n))"
              << "\n";

    // Space Complexity Row
    cout << setw(20) << "Space Complexity"
              << setw(15) << "O(n)"
              << setw(15) << "O(n)"
              << setw(15) << "O(n)"
              << "\n";
    
    cout<<endl;

    cout<<endl<<"The MergeSort begins here."<<endl;

    DivideEtImpera(arr, 0, arr.size() - 1);

    cout<<endl<<"And the MergeSort ends here."<<endl;

    cout <<endl<< "Sorted array: ";
    for(int i=0;i<n;i++)
        cout<<"On position "<<i<<" is the value "<<arr[i]<<"; arr["<<i<<"]="<<arr[i]<<endl;

    return 0;
}