#pragma once

class HMM;

class SimpleGenomicTest {
public:
	SimpleGenomicTest(Genome* g);

	HMM* genic, *intergenic;
};

