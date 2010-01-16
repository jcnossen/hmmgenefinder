#pragma once


template<class T>
class mvec
{
public:
	typedef int size_t;
	typedef T* ptr_t;
	typedef T* iterator;
	typedef const T* const_iterator;

	enum { min_size = 10 };

protected:
	// helpers for construction/destruction
	inline static void Constr (ptr_t ptr, const T& val) { new(ptr) T(val); }
	inline static void Destr (ptr_t ptr) { (ptr)->~T(); }

	// plain memory alloc/free, should not call destructors
	inline static ptr_t MAlloc (size_t n) { return (ptr_t)new char(n*sizeof(T)); }
	inline static void MFree (ptr_t p, size_t n) { char* c = (char*)p; delete[] c; } //(p,n*sizeof(T)); }

	inline static ptr_t Copy (ptr_t f, ptr_t l, ptr_t r) { for(;f!=l;++f,++r) *r=*f; return r; }
	inline static ptr_t UCopy (ptr_t f, ptr_t l, ptr_t r) { for(;f!=l;++f,++r) Constr(r,*f); return r; }
	inline static ptr_t Move (ptr_t f, ptr_t l, ptr_t r) { for(;f!=l;++f,++r) { Constr(r,*f); Destr(f); } return r; }
	inline static ptr_t RMove (ptr_t f, ptr_t l, ptr_t r) { for(r+=(l--)-f;l>=f;l--,r--) { Constr(r,*l); Destr(l); } return r; }
public:
	mvec() { First=Last=End=0; }
	mvec(const mvec<T>& c) { First=Last=End=0; copy_from(c); }
	~mvec() { erase(begin(),end()); MFree(First,End-First); }

	iterator begin() { return First; }
	iterator end() { return Last; }
	const_iterator begin() const { return First; }
	const_iterator end() const { return Last; }

	size_t size() const { return Last - First; }
	size_t capacity() const { return End - First; }
	void push_back(const T& v) { if(Last == End) realloc(); Constr(Last ++, v); }
	void pop_back() {assert(Last!=First);Destr(--Last);}
	void realloc() { size_t ns = First ? size()*2 : min_size; reserve(ns); }
	void clear() {erase(First,Last);}
	void reserve(size_t n) 
	{
		if (n < capacity ()) return;
		ptr_t New=MAlloc(n);
		if(size()) Move(First,Last,New);
		if(First) MFree(First,size());
		End=New+n;
		Last=New+size();
		First=New;
	}
	void erase(iterator F, iterator L)
	{
		for(iterator D=F;D!=L;++D) Destr(D);
		if(L<Last) Last=Move(L,Last,F);
		else Last=F;
	}
	void erase(iterator P) { erase(P,P+1); }
	void insert (iterator P, iterator F, iterator L)
	{
		if (End > Last - F + L) {
			RMove (P, Last, P - F + L);
			UCopy (F, L, P);
			Last += L-F;
		} else {
			ptr_t Ln, New = MAlloc (size() + (L - F));
			Ln = New;
			if (P > First) Ln = Move (First, P, Ln);
			Ln = UCopy (F, L, Ln);
			if (P < Last) Move (P, Last, Ln);
			MFree (First, size());
			End += New - First;
			Last = Ln;
			First = New;
		}
	}
	void resize(size_t n, const T& v = T()) 
	{
		if(size()<n) {
			reserve(n); for(;Last!=End;++Last) Constr(Last, v);
		} else erase (First+n,Last);
	}
	void copy_from(const mvec<T>& v)
	{
		if(&v==this) return;
		clear();reserve(v.size());
		Last=UCopy(v.First,v.Last,First);
	}
	bool empty() const { return Last==First; }
	T& back() { return *(Last-1); }
	const T& back() const { return *(Last-1); }
	T& front() { return *First; }
	const T& front() const { return *First; }

	// concatenation
	mvec operator&(const mvec& rhs)
	{
		mvec r(*this);
		r.add(rhs.begin(), rhs.end());
		return r;
	}

	// C++ style indexing
	T& operator[](size_t i) { return First[i]; }
	const T& operator[](size_t i) const { return First[i]; }
	// Matlab style indexing
	T& operator()(size_t i) { return First[i-1]; }
	const T& operator()(size_t i) const { return First[i-1]; }

	mvec<T>& operator=(const mvec<T>& v) { copy_from(v); return *this; }

	template<typename other_iterator> void add(other_iterator f, other_iterator l)
	{
		insert(end(), f, l);
	}

	template<typename R> mvec<R> cast() {
		mvec<R> r(size());
		for(int i=0;i<size();i++)
			r[i]=(R)First[i];
		return r;
	}

public:
	ptr_t First,
		Last, // Last points to the first unused entry (the one after the last entry)
		End;
};

template<typename T> T max(const mvec<T>& m) {
	if (m.empty())
		throw std::invalid_argument("max() called on empty mvec");

	T v=m.front();
	for(mvec<T>::const_iterator i=m.begin();i!=m.end();++i)
		v=std::max(v, *i); 
	return v;
}

template<typename T> T min(const mvec<T>& m) {
	if (m.empty())
		throw std::invalid_argument("min() called on empty mvec");

	T v=m.front();
	for(mvec<T>::const_iterator i=m.begin();i!=m.end();++i)
		v=std::min(v, *i); 
	return v;
}

