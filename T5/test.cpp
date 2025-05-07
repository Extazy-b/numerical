#include "./math.cpp"

using namespace std;

int main(){
    int N;
    cin >> N;
    vector <vector <size_t>> perms = permutations(N);
    for (vector <size_t> perm: perms)
    {
        for (size_t num: perm) {
            cout << num << ' ';
        }
        cout << endl;
    }
    return 0;
};