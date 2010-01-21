// dbdm4.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "Genome.h"

#include "math.h"


int main(int argc, char* argv[])
{
	try {
		Genome genome ("../data/AE005174.gd");

		genome.PrintInfo();

		mvec<Genome*> tt = genome.Split();
		Genome *train = tt(1), *test = tt(2);

		d_trace("Train genome: \n");
		train->PrintInfo();
		d_trace("Test genome: \n");
		test->PrintInfo();

		DeleteAll(tt);
	}
	catch (ContentException e) {
		d_trace("Exception: %s\n", e.what());
	}
	return 0;
}

