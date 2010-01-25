#pragma once

class HMM;
class Genome;

class SimpleGenomicTest {
public:
	SimpleGenomicTest(Genome* g);

	HMM* genic, *intergenic;
};

