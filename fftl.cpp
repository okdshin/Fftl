#ifdef FFTL_UNIT_TEST
#include "fftl.h"
//#include "Fftk/FastFourierTransform.h"
#include <iostream>

using namespace fftl;

template<class SignalArray>
auto print_signal_array(const std::string& mess, SignalArray& sa) -> void {
	std::cout << mess;
	using traits = signal_array_traits<SignalArray>;
	typename traits::iterator_type iter = traits::begin(sa);
	while(iter != traits::end(sa)){
		std::cout << *iter << " ";
		++iter;
	}
	std::cout << "\n";
}

int main(int argc, char* argv[])
{
	int n = 128;
	{
		fast_basic_transform<std::complex<double>, 
			std::vector<std::complex<double>>> fft(n);

		std::vector<std::complex<double>> src;
		for(unsigned int i = 0; i < n; ++i){
			src.push_back(std::complex<double>(sin(M_PI*(static_cast<double>(i)/n)), 0));
		}
		print_signal_array("src: ", src);

		auto transformed = fft.transform(src);
		print_signal_array("tra: ", transformed);

		auto inversed = fft.inverse_transform(transformed);
		print_signal_array("inv: ", inversed);
	}    

	/*
	std::cout << "\n";

	{
		using namespace fftk;
		FastFourierTransform fft((SignalLength(n)));
		
		std::vector<std::complex<double>> src;
		for(unsigned int i = 0; i < n; ++i){
			src.push_back(std::complex<double>(sin(M_PI*(static_cast<double>(i)/n)), 0));
		}

		auto transformed = fft.Transform(src);
		//print_signal_array("tra: ", transformed);
		
		auto inversed = fft.InverseTransform(transformed);
		print_signal_array("inv: ", inversed);

		
	}
	*/
	return 0;
}

 #endif
