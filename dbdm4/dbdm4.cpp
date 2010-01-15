// dbdm4.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "Genome.h"

#include "math.h"

int main(int argc, char* argv[])
{

	try {
		Genome genome ("..\\data\\AE005174.gd");

		genome.PrintInfo();
	}
	catch (ContentException e) {
		d_trace("Exception: %s\n", e.what());
	}
	return 0;
}

