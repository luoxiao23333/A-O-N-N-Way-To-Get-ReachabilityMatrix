#include<iostream>
#include<vector>
#include<stack>
#include<algorithm>
#include<set>

using namespace std;

/*
定义 N为图的节点数，E为边的数量(除去数组初始化)
步骤						时间复杂度						空间复杂度
邻接矩阵					O(E)							O(N*N)
关联矩阵					O(E)							O(E*N)
强分图						O(E)							O(E)
可达矩阵					O(N*N)							O(N*N)
*/

//矩阵基类
class Matrix {
protected:

	int**Data;

	int ColLength,
		RawLength;

	constexpr void SetValue(const int&col, const int&raw, const int&&val) {
		Data[col][raw] = val;
	}


public:

	const int&GetColLength() { return this->ColLength; }

	const int&GetRawLength() { return this->RawLength; }

	Matrix(const int&Col, const int&Raw) :ColLength(Col), RawLength(Raw) {
		Data = new int*[Col + 1];
		for (int col = 1; col <= Col; ++col) {
			Data[col] = new int[Raw + 1];
			for (int raw = 1; raw <= Raw; ++raw) {
				Data[col][raw] = 0;
			}
		}
	}

	virtual void Output() {
		cout << '\t';
		for (int i = 1; i <= ColLength; ++i) {
			cout << 'v' << i << '\t';
		}
		cout << endl;
		for (int col = 1; col <= ColLength; ++col) {
			cout << 'v' << col << '\t';
			for (int raw = 1; raw <= RawLength; ++raw) {
				cout << Data[col][raw] << '\t';
			}
			cout << endl;
		}
	}

	~Matrix() {
		for (int i = 1; i <= ColLength; ++i) {
			delete[]Data[i];
		}
		delete[]Data;
	}
};

//邻接矩阵
class AdjointMatrix :public Matrix {

	friend class StrongPartiteGraph;
	friend class ReachabilityMatrix;

public:

	AdjointMatrix(const int&PointNumber) :Matrix(PointNumber, PointNumber) {  }

	void AddEdge(const int&from, const int&to) {
		this->SetValue(from, to, 1);
	}

	constexpr bool IsConnected(const int&from, const int&to) {
		return Data[from][to];
	}

	void Output()override {
		cout << '\t';
		for (int i = 1; i <= ColLength; ++i) {
			cout << 'v' << i << '\t';
		}
		cout << endl;
		for (int col = 1; col <= ColLength; ++col) {
			cout << 'v' << col << '\t';
			for (int raw = 1; raw <= RawLength; ++raw) {
				cout << Data[col][raw] << '\t';
			}
			cout << endl;
		}
	}

	AdjointMatrix*Clone() {
		AdjointMatrix*AM = new AdjointMatrix(this->ColLength);
		size_t&&Size = sizeof(int)*(ColLength + 1);
		for (int Index = 1; Index <= this->ColLength; ++Index) {
			memcpy(AM->Data[Index], this->Data[Index], Size);
		}
		return AM;
	}

};

//关联矩阵
class RelatedMatrix :public Matrix {
public:

	RelatedMatrix(const int&PointNumber, const int&EdgeNumber) :Matrix(PointNumber, EdgeNumber) {  }

	void AddEdge(const int&EdgeIndex, const int&from, const int&to) {
		this->SetValue(from, EdgeIndex, 1);
		this->SetValue(to, EdgeIndex, -1);
	}

	void Output()override {
		cout << '\t';
		for (int i = 1; i <= RawLength; ++i) {
			cout << 'e' << i << '\t';
		}
		cout << endl;
		for (int col = 1; col <= ColLength; ++col) {
			cout << 'v' << col << '\t';
			for (int raw = 1; raw <= RawLength; ++raw) {
				cout << Data[col][raw] << '\t';
			}
			cout << endl;
		}
	}

};

//图
class Graph {

	vector<int>*G;

	struct Edge {
		const int from,to;
		Edge(const int&from, const int&to) :from(from), to(to) {  }
	};
	vector<Edge> Edges;

	int PointNumber;

	friend class StrongPartiteGraph;

public:

	Graph(const int&N) :PointNumber(N) { 
		G = new vector<int>[PointNumber + 1];
	}

	void AddEdge(const int&from, const int&to) {
		G[from].push_back(to);
		Edges.emplace_back(from, to);
	}

	AdjointMatrix* GetAdjointMatrix() {
		AdjointMatrix*M = new AdjointMatrix(PointNumber);
		for (int from = 1; from <= PointNumber; ++from) {
			for (const auto&to : G[from]) {
				M->AddEdge(from, to);
			}
		}
		return M;
	}

	RelatedMatrix* GetRelatedMatrix() {
		int&&EdgeNumber = 0;
		for (int from = 1; from <= PointNumber; ++from) {
			EdgeNumber += G[from].size();
		}
		RelatedMatrix*M = new RelatedMatrix(PointNumber, EdgeNumber);
		for (int Index = 0; Index < static_cast<int>(Edges.size()); ++Index) {
			M->AddEdge(Index + 1, Edges[Index].from, Edges[Index].to);
		}
		return M;
	}

	~Graph() {
		delete[] G;
	}

};

//并查集，用于算法优化
class DisjointSet {

	int MaxPoint;

	int *Father;

public:

	DisjointSet(int MaxPoint) {
		this->Father = new int[MaxPoint + 1];
		for (int i = 1; i <= MaxPoint; ++i) {
			this->Father[i] = i;
		}
	}

	int Find(int x) {
		return this->Father[x] == x ? x : (this->Father[x] = Find(this->Father[x]));
	}

	void Merge(int x, int y) {
		this->Father[this->Find(x)] = this->Find(y);
	}

	bool isSameSet(int x, int y) {
		return Find(x) == Find(y);
	}

	~DisjointSet() {
		delete[]this->Father;
	}

};

//强分图,采用Gabow算法
class StrongPartiteGraph {

	vector<int> *G;

	AdjointMatrix*StrongPartiteMatrix;

	vector<int>*PointSet;

	vector<int>*StrongPartiteList;

	const int&PointNumber;

	int *Belongs,
		Color,
		*TimeTag,
		*Stack1,
		*Stack2;

	bool *Visable;

	void InitAdjointMatrix() {

		int&&CurrentSet = 1;

		this->PointSet = new vector<int>[Color + 1];
		this->StrongPartiteList = new vector<int>[Color + 1];

		this->StrongPartiteMatrix = new AdjointMatrix(Color);

		for (int Index = 1; Index <= PointNumber; ++Index) {
			this->PointSet[this->Belongs[Index]].push_back(Index);
		}

		DisjointSet DS(PointNumber);

		for (int Index = 1; Index <= Color; ++Index) {

			if (this->PointSet[Index].size() == 1)continue;

			auto it = this->PointSet[Index].begin();
			const int&temp = *(++it);
			for (; it != this->PointSet[Index].end(); ++it) {
				DS.Merge(*it, temp);
			}

		}

		for (int Index = 1; Index <= Color; ++Index) {
			for (const auto&Point : this->PointSet[Index]) {
				for (const auto&to : this->G[Point]) {
					if (!DS.isSameSet(Point, to)) {
						this->StrongPartiteMatrix->AddEdge(Index, this->Belongs[to]);
						if (Index != this->Belongs[to]) {
							this->StrongPartiteList[Index].push_back(this->Belongs[to]);
						}
					}
				}
			}
		}

	}

	void Visit(int cur, int&sig, int&scc_num){

		TimeTag[cur] = ++sig;

		Stack1[++Stack1[0]] = cur;
		Stack2[++Stack2[0]] = cur;

		for (const auto&to : this->G[cur]){
			if (!TimeTag[to]){
				Visit(to, sig, scc_num);
			}
			else if (!Belongs[to]){
				while (TimeTag[Stack2[Stack2[0]]] > TimeTag[to]) {
					--Stack2[0];
				}
			}
		}

		if (Stack2[Stack2[0]] == cur){

			--Stack2[0]; 
			++scc_num;

			do{
				Belongs[Stack1[Stack1[0]]] = scc_num;
			} while (Stack1[Stack1[0]--] != cur);

		}
	}

	void GabowCompute(){
		int sig;
		this->TimeTag = new int[PointNumber + 1];
		this->Stack1 = new int[PointNumber + 1];
		this->Stack2 = new int[PointNumber + 1];

		memset(this->Belongs, 0, sizeof(int)*(PointNumber + 1));
		memset(TimeTag, 0, sizeof(int)*(PointNumber + 1));

		sig = 0; 
		Color = 0;
		Stack1[0] = 0; 
		Stack2[0] = 0;

		for (int Point = 1; Point <= PointNumber; ++Point){
			if (!TimeTag[Point])
				Visit(Point, sig, Color);
		}

		delete[]this->TimeTag;
		delete[]this->Stack1;
		delete[]this->Stack2;

		this->InitAdjointMatrix();

	}

public:
	StrongPartiteGraph(const Graph&G) :G(G.G), PointNumber(G.PointNumber) {

		this->Belongs = new int[PointNumber + 1];
		this->Visable = new bool[PointNumber + 1];

		memset(this->Belongs, 0x0, sizeof(int)*(PointNumber + 1));
		memset(this->Visable, 0x0, sizeof(bool)*(PointNumber + 1));

		this->GabowCompute();

	}

	AdjointMatrix*GetReachabilityMatrix() {

		AdjointMatrix*AM = new AdjointMatrix(PointNumber);
		int cnt = 0;

		auto Translate = [&AM,this](const int&FromSetIndex, const int&ToSetIndex)->void {
			for (const auto&from : this->PointSet[FromSetIndex]) {
				for (const auto&to : this->PointSet[ToSetIndex]) {
					AM->AddEdge(from, to); 
				}
			}
		};

		bool**Vis = new bool*[PointNumber + 1];

		for (int Index = 1; Index <= PointNumber; ++Index) {
			Vis[Index] = new bool[PointNumber + 1];
			memset(Vis[Index], 0x0, sizeof(bool)*(PointNumber + 1));
		}

		stack<int>Stack;
		for (int Point = 1; Point <= Color; ++Point) {
			Stack.push(Point);
			vector<int>Froms;
			while (!Stack.empty()) {
				int from = Stack.top();
				Stack.pop();
				Froms.push_back(from);
				for (const auto&Precursor : Froms) {
					if (Vis[Precursor][from])break;
					Vis[Precursor][from] = true;
					Translate(Precursor, from);
				}
				for (const auto&to : this->StrongPartiteList[from]) {
					Stack.push(to);
				}
			}
		}

		for (int Index = 1; Index <= PointNumber; ++Index) {
			delete[]Vis[Index];
		}
		delete[]Vis;

		return AM;
	}

	void Output() {

		cout << "强分图共有" << Color << "个强连通分量" << endl
			<< "强分图为" << endl;
		this->StrongPartiteMatrix->Output();
		cout << "其中：" << endl;
		for (int Index = 1; Index <= Color; ++Index) {
			cout << "v" << Index << " = {";
			for (const auto&Point : this->PointSet[Index]) {
				cout << Point << ',';
			}
			cout << '\b' << '}' << endl;
		}
	}

	~StrongPartiteGraph() {
		delete[]this->Belongs;
		delete[]this->PointSet;
		delete this->StrongPartiteMatrix;
		delete[]this->Visable;
	}

};

int main() {

	int N;
	cout << "请输入点的数量N（编号从1到N）:" << endl;
	cin >> N;

	Graph G(N);

	int M;
	cout << "请输入有向边的数量M:" << endl;
	cin >> M;

	cout << "依次输入有向边:" << endl
		<< "from\tto" << endl;
	while (M--) {
		int from, to;
		cin >> from >> to;
		G.AddEdge(from, to);
	}

	Matrix* matrix;
	cout << "邻接矩阵：" << endl;
	matrix = G.GetAdjointMatrix();
	matrix->Output();
	delete matrix;

	cout << "关联矩阵：" << endl;
	matrix = G.GetRelatedMatrix();
	matrix->Output();
	delete matrix;

	cout << "强分图：" << endl;
	StrongPartiteGraph SPG(G);
	SPG.Output();

	cout << "可达矩阵：" << endl;
	matrix = SPG.GetReachabilityMatrix();
	matrix->Output();
	delete matrix;

	system("pause");

	return 0;

}
