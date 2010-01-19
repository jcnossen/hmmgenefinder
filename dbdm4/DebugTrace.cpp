#include "stdafx.h"
#include <cstdarg>

#if defined(_MSC_VER) && defined(_DEBUG)
#define DEBUGTRACE
#endif

#ifdef DEBUGTRACE
	#include <windows.h>
#endif

static std::string logFile;

void d_trace(const char *fmt, ...)
{
	static char buf[256];

	va_list ap;
	va_start(ap,fmt);
	VSNPRINTF (buf,256,fmt,ap);
	va_end(ap);

	d_puts (buf);
}

/*
d_puts - writes a string of any length to the debug output
	(as opposed to d_trace, which is limited by the buffer)
*/
void d_puts (const char *buf)
{
#ifdef DEBUGTRACE
	OutputDebugStringA(buf);
#endif

	fputs(buf, stderr);
	
	if (!logFile.empty()) {
		// assuming app can crash at any time, its better to open and close then to keep the file open.
		FILE *f = fopen(logFile.c_str(), "a");
		fputs(buf, f);
		fclose(f);
	}
}



void d_assert(char *str, char *file, int line)
{
	std::string msg = SPrintf("Assertion \"%s\" failed at line %d in \'%s\'\n", str, line, file);
	d_trace (msg.c_str());

	throw std::runtime_error(msg.c_str());
}

void d_setlogfile( const char *f )
{
	logFile = f;
}
