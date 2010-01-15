#pragma once


// Matlab style vector
template<typename T>
class mvec {
	T* d; // elements
	int n;
public:
	mvec() : n(0),d(0) { }
	mvec(int size) {
		n=size;
		d=new T[size];
	}
	~mvec() {
		delete d;
	}

	T& operator[](int i) { return d[i-1]; }
	const T& operator[](int i) const { return d[i-1]; }
	int size() { return n; }

	static mvec range(T a, T b) {
		return range(a, 1, b);
	}
	static mvec range(T a, T step, T b) {
		int n = 1+(b-a)/step;
		mvec r(n);
		for(int i=0;i<n;i++)
			r.d[i] = a+i*step;
		return r;
	}
};

