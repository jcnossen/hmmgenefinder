/* 
Matlab style vector class
Copyright 2010, Jelmer Cnossen


This class tries to support matlab style vectors by

- allowing +-/* operators on combinations of vector,vector and vector,scalar
- using the () operators for indexing with first index 1
- implementing max, min, mean, sum functions that operate on the mvec instances

However since there is no .* operator in C++, inproducts require regular functions and can't be done with *
Also for similar reasons, concatenation is done with & operator instead of [a b]
*/
#pragma once

template<typename T> class type_traits {};
template<> class type_traits<float> {
public:
	static float zero() { return 0.0f; }
};
template<> class type_traits<int> {
public:
	static int zero() { return 0; }
};
template<> class type_traits<double> {
public:
	static double zero() { return 0.0; }
};
template<typename T> class type_traits<T*> {
public:
	static T* zero() { return NULL; }
	typedef T value_type;
};

template<typename T>
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
	inline static ptr_t MAlloc (size_t n) { return (ptr_t)new char[n*sizeof(T)]; }
	inline static void MFree (ptr_t p, size_t n) { char* c = (char*)p; delete[] c; } //(p,n*sizeof(T)); }

	// use operator=
	inline static ptr_t Copy (ptr_t f, ptr_t l, ptr_t r) { for(;f!=l;++f,++r) *r=*f; return r; }
	// use copy contructor (for uninitialized memory)
	template<typename TSrc> inline static ptr_t UCopy (TSrc* f, TSrc* l, ptr_t r) { for(;f!=l;++f,++r) Constr(r,*f); return r; }
	inline static ptr_t Move (ptr_t f, ptr_t l, ptr_t r) { for(;f!=l;++f,++r) { Constr(r,*f); Destr(f); } return r; }
	inline static ptr_t RMove (ptr_t f, ptr_t l, ptr_t r) { for(r+=(l--)-f;l>=f;l--,r--) { Constr(r,*l); Destr(l); } return r; }
public:
	mvec() { First=Last=End=0; }
	mvec(const mvec<T>& c) { First=Last=End=0; set(c); }
	~mvec() { erase(begin(),end()); MFree(First,End-First); First=End=Last=0; }

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
	template<typename SrcIterator>
	void insert (iterator P, SrcIterator F, SrcIterator L)
	{
		if (End - Last > L - F) { // room left for new items
			RMove (P, Last, P + (L - F)); // move [P, Last) to P + (L-F). IAW make room for new items
			UCopy (F, L, P);
			Last += L-F;
		} else {
			int n = size() + (L-F);
			ptr_t Ln, New = MAlloc (n);
			Ln = New;
			if (P > First) Ln = Move (First, P, Ln);
			Ln = UCopy (F, L, Ln);
			if (P < Last) Ln = Move (P, Last, Ln);
			if (First) MFree (First, size());
			Last = End = New + n; // Ln = End+n
			First = New;
		}
	}
	void resize(size_t n, const T& v = T()) 
	{
		if(size()<n) {
			reserve(n); for(;Last!=End;++Last) Constr(Last, v);
		} else erase (First+n,Last);
	}
	template<typename TOther>
	void set(const mvec<TOther>& v)
	{
		if((void*)&v==(void*)this) 
			return;
		clear();reserve(v.size());
		Last=UCopy(v.First,v.Last,First);
	}
	bool empty() const { return Last==First; }
	T& back() { return *(Last-1); }
	const T& back() const { return *(Last-1); }
	T& front() { return *First; }
	const T& front() const { return *First; }

	// concatenation
	template<typename TSrcContainer>
	mvec operator&(const TSrcContainer& rhs) const
	{
		mvec r(*this);
		r.add(rhs.begin(), rhs.end());
		return r;
	}
	template<typename TSrcContainer>
	mvec& operator&=(const TSrcContainer& rhs) {
		add(rhs.begin(), rhs.end());
		return *this;
	}

	// C++ style indexing
	T& operator[](size_t i) { return First[i]; }
	const T& operator[](size_t i) const { return First[i]; }
	// Matlab style indexing
	T& operator()(size_t i) { return First[i-1]; }
	const T& operator()(size_t i) const { return First[i-1]; }

	// For reasons beyond my C++ skills, we need 2 =operators here
	mvec<T>& operator=(const mvec<T>& v) { set(v); return *this; }
	template<typename B> mvec<T>& operator=(const mvec<B>& v) { set(v); return *this; }

	template<typename other_iterator> void add(other_iterator f, other_iterator l)
	{
		insert(end(), f, l);
	}

	template<typename R> mvec<R> cast() const {
		mvec<R> r(size());
		for(int i=0;i<size();i++)
			r[i]=(R)First[i];
		return r;
	}

	// Assumes object container
	mvec clone() const {
		mvec r;
		r.reserve(size());
		for(iterator i=First;i!=Last;++i)
			r.push_back(new typename type_traits<T>::value_type(**i));
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
	for(typename mvec<T>::const_iterator i=m.begin();i!=m.end();++i)
		v=std::max(v, *i); 
	return v;
}

template<typename T> T min(const mvec<T>& m) {
	if (m.empty())
		throw std::invalid_argument("min() called on empty mvec");

	T v=m.front();
	for(typename mvec<T>::const_iterator i=m.begin();i!=m.end();++i)
		v=std::min(v, *i); 
	return v;
}
template<typename T> T sum(const mvec<T>& m) {
	if (m.empty())
		throw std::invalid_argument("sum() called on empty mvec");

	T v=type_traits<T>::zero;
	for(typename mvec<T>::const_iterator i=m.begin();i!=m.end();++i)
		v+=*i;
	return v;
}
template<typename T> float mean(const mvec<T>& m) {
	return sum(m) / (float)m.size();
}

namespace ops {
	template<typename A, typename B> struct add { static A apply(A a, B b) { return a+b; } 	};
	template<typename A, typename B> struct sub { static A apply(A a, B b) { return a-b; } 	};
	template<typename A, typename B> struct mul { static A apply(A a, B b) { return a*b; } 	};
	template<typename A, typename B> struct div { static A apply(A a, B b) { return a/b; } 	};
	// reverse order apply
	template<typename A, typename B, typename Op> struct rev { static A apply(A a, B b) { return Op::apply(b, a); } };
};

template<typename TA, typename TB, typename TOperatorType> 
class operator_helper {
public:
	static mvec<TA> apply(const mvec<TA>& a, TB b) {
		mvec<TA> r; r.reserve(a.size());
		for(typename mvec<TA>::const_iterator p=a.begin();p!=a.end();++p)
			r.push_back( TOperatorType::apply((*p), b) );
		return r;
	}
};
template<typename TA, typename TVec, typename TOperatorType>
class operator_helper<TA, const mvec<TVec>&, TOperatorType > { // specialization for container
public:
	static mvec<TA> apply(const mvec<TA>& a, const mvec<TVec>& b) {
		mvec<TA> r; r.reserve(a.size());
		if (b.size() != a.size()) 
			throw std::invalid_argument("mvec sizes do not match for elementwise operation");
		for (int x=0;x<a.size();x++)
			r.push_back(TOperatorType::apply(a[x], b[x]));
		return r;
	}
};

#define OPERATOR_IMPL(SYM, OP) \
template<typename T> static mvec<T> operator SYM(const mvec<T>& c, const mvec<T>& v) { \
	return operator_helper<T, const mvec<T>&, ops::OP<T,T> >::apply(c,v); \
}  \
template<typename T, typename B> static mvec<T> operator SYM(const mvec<T>& c, B v) { \
	return operator_helper<T, B, ops::OP<T, B> >::apply(c, v);  \
} \
template<typename T, typename B> static mvec<T> operator SYM(B v, const mvec<T>& c) { \
	typedef ops::rev<T, B, ops::OP<B,T> > operator_t;  \
	return operator_helper<T, B, operator_t>(c,v); \
}

OPERATOR_IMPL(+, add)
OPERATOR_IMPL(-, sub)
OPERATOR_IMPL(*, mul)
OPERATOR_IMPL(/, div)

#undef OPERATOR_IMPL


typedef mvec<float> mvecf;
typedef mvec<int> mveci;
typedef mvec<double> mvecd;
