/* 
Matlab style vector class
Copyright 2010, Jelmer Cnossen


This class tries to support matlab style vectors by

- allowing +-/* operators on combinations of vector,vector and vector,scalar
- using the () operators for indexing with first index 1
- implementing max, min, mean, sum functions that operate on the mvec instances
- members() and members_ptr() function, to replace matlabs [array.member]

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

//default cloning: just use copy constructor.
// use template specialization to create other methods of cloning
template<typename T> class cloning_device {
public:
	static T* clone(T* src) { return new T(*src); }
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
	mvec(size_t s, const T& v = T()) { First=Last=End=0; resize(s, v); }
	mvec(const mvec<T>& c) { First=Last=End=0; set(c); }
	mvec(T* first, size_t count) { First=Last=End=0; insert(0, first, first+count); }
	template<typename TSrc>
	mvec(TSrc* begin_, TSrc* end_) { First=Last=End=0; insert(0, begin_, end_); }
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
		if (!n || n < capacity ()) return;
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

	template<typename B>
	mvec<T> operator-() {
		mvec<T> r; r.reserve(size());
		for(int i=0;i<size();i++)
			r[i]=-First[i];
		return r;
	}

	mvec<int> operator!() {
		mvec<T> r; r.reserve(size());
		for(int i=0;i<size();i++)
			r[i]=First[i] ? 0 : 1;
		return r;
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
			r.push_back(cloning_device<type_traits<T>::value_type>::clone(*i));
		return r;
	}

	// overloaded -> to assign to set of pointers
	template<typename V>
	void ptr_assign(const mvec<V>& vals) {
		assert(vals.size()==size());
		for(int i=0;i<size();i++)
			*(First[i])=vals[i];
	}
	template<typename V>
	void ptr_assign(const V& val) {
		for(int i=0;i<size();i++)
			*(First[i])=val;
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

	T v=type_traits<T>::zero();
	for(typename mvec<T>::const_iterator i=m.begin();i!=m.end();++i)
		v+=*i;
	return v;
}
template<typename T> float mean(const mvec<T>& m) {
	return sum(m) / (float)m.size();
}

template<typename T, typename T2, typename TMember> mvec<TMember> members(typename const mvec<T*>& src, TMember T2::*member) {
	mvec<TMember> r; r.reserve(src.size());
	for(int i=0;i<src.size();i++)
		r.push_back(src[i]->*member);
	return r;
}

template<typename T, typename T2, typename TMember> mvec<TMember> members(typename const mvec<T>& src, TMember T2::*member) {
	mvec<TMember> r; r.reserve(src.size());
	for(int i=0;i<src.size();i++)
		r.push_back(src[i].*member);
	return r;
}

template<typename T, typename T2, typename TMember> mvec<TMember*> members_ptr(typename mvec<T*>& src, TMember T2::*member) {
	mvec<TMember*> r; r.reserve(src.size());
	for(int i=0;i<src.size();i++)
		r.push_back(&( (src[i])->*member ));
	return r;
}

template<typename T, typename T2, typename TMember> mvec<TMember*> members_ptr(typename mvec<T>& src, TMember T2::*member) {
	mvec<TMember*> r; r.reserve(src.size());
	for(int i=0;i<src.size();i++)
		r.push_back(&src[i].*member);
	return r;
}

template<typename T>
typename mvec<T>& operator|=(mvec<T>& ptrs, const mvec< typename type_traits<T>::value_type >& values ) {
	ptrs.ptr_assign(values);
	return ptrs;
}
template<typename T, typename B>
typename mvec<T>& operator|=(mvec<T>& ptrs, const B& val) {
	ptrs.ptr_assign(val);
	return ptrs;
}


namespace mvec_ops {
	template<typename A, typename B> struct add { 
		typedef A return_type;
		static A apply(A a, B b) { return a+b; } 	
	};
	template<typename A, typename B> struct sub { 
		typedef A return_type;
		static A apply(A a, B b) { return a-b; } 	
	};
	template<typename A, typename B> struct mul { 
		typedef A return_type;
		static A apply(A a, B b) { return a*b; } 	
	};
	template<typename A, typename B> struct div { 
		typedef A return_type;
		static A apply(A a, B b) { return a/b; } 	
	};
	// reverse order apply
	template<typename A, typename B, typename Op> struct rev { 
		typedef typename Op::return_type return_type;
		static A apply(A a, B b) { return Op::apply(b, a); } 
	};
	template<typename A, typename B> struct eq {
		typedef int return_type;
		static int apply(A a, B b) { return a==b ? 1 : 0; }
	};
	template<typename A, typename B> struct neq {
		typedef int return_type;
		static int apply(A a, B b) { return a!=b ? 1 : 0; }
	};
	template<typename A, typename B> struct ge {
		typedef int return_type;
		static int apply(A a, B b) { return a>=b ? 1 : 0; }
	};
	template<typename A, typename B> struct le {
		typedef int return_type;
		static int apply(A a, B b) { return a<=b ? 1 : 0; }
	};
	template<typename A, typename B> struct l_or {
		typedef int return_type;
		static int apply(A a, B b) { return a||b ? 1 : 0; }
	};
	template<typename A, typename B> struct l_and {
		typedef int return_type;
		static int apply(A a, B b) { return a&&b ? 1 : 0; }
	};
};

template<typename TA, typename TB, typename TOperatorType> 
class operator_helper {
public:
	typedef typename TOperatorType::return_type RT;
	static mvec<RT> apply(const mvec<TA>& a, TB b) {
		mvec<RT> r; r.reserve(a.size());
		for(typename mvec<TA>::const_iterator p=a.begin();p!=a.end();++p)
			r.push_back( TOperatorType::apply((*p), b) );
		return r;
	}
};

template<typename TA, typename TVec, typename TOperatorType>
class operator_helper<TA, const mvec<TVec>&, TOperatorType > { // specialization for container
public:
	typedef typename TOperatorType::return_type RT;
	static mvec<RT> apply(const mvec<TA>& a, const mvec<TVec>& b) {
		mvec<RT> r; r.reserve(a.size());
		if (b.size() != a.size()) 
			throw std::invalid_argument("mvec sizes do not match for elementwise operation");
		for (int x=0;x<a.size();x++)
			r.push_back(TOperatorType::apply(a[x], b[x]));
		return r;
	}
};


#define OPERATOR_IMPL(SYM, OP) \
template<typename T> inline typename mvec< typename mvec_ops::OP<T,T>::return_type > operator SYM(const mvec<T>& c, const mvec<T>& v) {	\
	return operator_helper<T, const mvec<T>&, mvec_ops::OP<T,T> >::apply(c,v);	\
}	\
template<typename T, typename B> inline typename mvec< typename mvec_ops::OP<T,B>::return_type > operator SYM(const mvec<T>& c, B v) {	\
	return operator_helper<T, B, mvec_ops::OP<T, B> >::apply(c, v);	\
}	\
template<typename T, typename B> inline typename mvec< typename mvec_ops::OP<T,B>::return_type> operator SYM(B v, const mvec<T>& c) {	\
	typedef mvec_ops::rev<T, B, mvec_ops::OP<B,T> > operator_t;	\
	return operator_helper<T, B, operator_t>(c,v);	\
}

OPERATOR_IMPL(+, add)
OPERATOR_IMPL(-, sub)
OPERATOR_IMPL(*, mul)
OPERATOR_IMPL(/, div)
OPERATOR_IMPL(==, eq)
OPERATOR_IMPL(!=, neq)
OPERATOR_IMPL(>=, ge)
OPERATOR_IMPL(<=, le)
OPERATOR_IMPL(||, l_or)
OPERATOR_IMPL(&&, l_and)

#undef OPERATOR_IMPL

#define OPERATOR_IMPL(SYM) \
template<typename A, typename B> inline mvec<A>& operator SYM(mvec<A>& a, const mvec<B>& b) { \
	if (b.size() != a.size())		\
		throw std::invalid_argument("mvec sizes do not match for elementwise operation"); \
	for (int i=0;i<a.size();i++) \
		a[i] SYM b[i]; \
	return a; \
} \
template<typename A, typename B> \
inline mvec<A>& operator SYM(mvec<A>& a, const B& b) { \
	for (int i=0;i<a.size();i++) \
		a[i] SYM b; \
	return a; \
}

OPERATOR_IMPL(+=)
OPERATOR_IMPL(-=)
OPERATOR_IMPL(*=)
OPERATOR_IMPL(/=)

#undef OPERATOR_IMPL


typedef mvec<float> mvecf;
typedef mvec<int> mveci;
typedef mvec<double> mvecd;
