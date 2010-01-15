#ifndef JC_CONFIG_PARSER_H
#define JC_CONFIG_PARSER_H

#include <string>
#include <list>
#include <cstdio>
#include <vector>

class CfgValue;
class CfgListElem;
class CfgList;
class CfgNumeric;
class CfgLiteral;


struct InputBuffer;


class CfgWriter
{
public:
	CfgWriter(const char *name);
	~CfgWriter();
	bool HasFailed();

	CfgWriter& operator<<(const std::string& s);
	CfgWriter& operator<<(char c);
	CfgWriter& operator<<(const char *str);

	void MakeIndent (char c);

	void IncIndent();
	void DecIndent();
protected:
	int indentLevel;
	FILE *out;
};

class CfgValue
{
public:
	enum Type {
		List, Numeric, Literal
	};

	CfgValue() {}
	virtual ~CfgValue() {}

	virtual void Parse (InputBuffer& buf) = 0;
	virtual void dbgPrint (int depth);
	virtual void Write (CfgWriter& w) = 0;
	virtual Type getType() = 0;

	std::string getFilePath() { return cfgFilePath; }

	static char Lookahead (InputBuffer& buf);
	static CfgValue* ParseValue(InputBuffer& buf);
	static CfgList* LoadNestedFile(InputBuffer& buf);
	static CfgList* LoadFile(std::string name);

protected:

	std::string cfgFilePath;
};

class CfgListElem
{
public:
	CfgListElem(CfgValue* v, const char* n=0) { value=v; if(n) name=n; }
	~CfgListElem() { if (value) delete value; }
	void Parse (InputBuffer& buf);
	void Write (CfgWriter& w);

	std::string name;
	CfgValue *value;
};

// list of cfg nodes
class CfgList : public CfgValue
{
public:
	~CfgList();

	void Parse (InputBuffer& buf, bool root);
	void Parse (InputBuffer& buf) { return Parse (buf, false); }
	void Write (CfgWriter& w) { Write(w,false); }
	void Write (CfgWriter& w, bool root);
	Type getType () { return CfgValue::List; }

	void dbgPrint (int depth);

	CfgValue* GetValue (const char *name);
	double GetNumeric(const char *name, double def=0.0);
	const char* GetLiteral(const char *name, const char *def=0);
	int GetInt (const char *name, int def=0) { return (int) GetNumeric (name, def); }
	CfgList* GetList(const char *name);

	CfgValue* Index(int i);

	static CfgList* EmptyList;

	void Add (const char *name, const char *val);
	void Add (const char *name, const std::string& s) { Add(name, s.c_str()); }
	void Add (const char *name, double val);
	void Add (const char *name,CfgValue *val);

	void Add (const std::string& s) { Add(0, s); }
	void Add (const char *val) { Add(0, val); }
	void Add (double val) { Add(0, val); }
	void Add (CfgValue *val) { Add(0, val); }

	typedef std::vector<CfgListElem*> ChildList;
	ChildList childs;
	typedef ChildList::iterator iterator;

	iterator begin() { return childs.begin(); }
	iterator end() { return childs.end(); }

	std::string& Load(const char *name, std::string& val) { const char *r = GetLiteral(name); if (r) val=r; return val; }
	CfgValue*& Load(const char *name, CfgValue*& val) { val=GetValue(name); return val; }
	bool& Load(const char* name, bool& val) { val=GetNumeric(name)!=0.0f; return val; }
	template<typename T> T& LoadNVal(const char *name, T& val) { val=(T)GetNumeric(name); return val; }
};

class CfgLiteral : public CfgValue  {
public:
	CfgLiteral () { }
	void Parse (InputBuffer& buf);
	void Write (CfgWriter& w);
	void dbgPrint (int depth);

	Type getType() { return CfgValue::Literal; }

	std::string value;
};

class CfgNumeric : public CfgValue {
public:
	CfgNumeric () { value=0.0; }
	void Parse (InputBuffer& buf);
	void Write (CfgWriter& w);
	void dbgPrint (int depth);
	Type getType() { return CfgValue::Numeric; }

	double value;
};


#endif
