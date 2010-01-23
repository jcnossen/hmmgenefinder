
#pragma once


 // same as matlab: A => 1 , C => 2, G => 3,  T(U) => 4.
int nt2int(char v);
mvec<int> nt2int(std::string v);

std::string seqcomplement(std::string dna);
std::string segrcomplement(std::string dna);
