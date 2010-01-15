
#pragma once


struct InputBuffer
{
	InputBuffer () : pos(0), line(1), data(0), len(0), filename(0) {}
	char operator*() const { return data[pos]; }
	bool end() const {  return pos == len; }
	InputBuffer& operator++() { next(); return *this; }
	InputBuffer operator++(int) { InputBuffer t=*this; next(); return t; }
	char operator[](int i)const { return data[pos+i]; }
	void showLocation() const;
	bool compareIdent (const char *str) const;
	char get() const { return data[pos]; }
	char get(int i) const { return data[pos+i]; }
	void next() { if (pos<len) pos++; }
	bool skipWhitespace();
	std::string location() const;
	void skipKeyword(const char *s);
	void expecting(const char *s); // show an error message saying that 's' was expected
	std::string parseIdent();

	int pos;
	int line;
	char *data;
	int len;
	const char *filename;
};


std::string ReadZStr(FILE *f);
void WriteZStr(FILE *f, const std::string& s);
std::string SPrintf(const char *fmt, ...);
std::string GetDirPath(std::string file); // strip filename
std::string ReadTextFile(std::string file);
