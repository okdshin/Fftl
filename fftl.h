#ifndef FFTL_ROTOR_TRAITS_H
#define FFTL_ROTOR_TRAITS_H

#include <iostream>
#include <stdexcept>
#include <complex>
#include <vector>

namespace fftl {

	template<class Rotor>
	struct rotor_traits {
		using value_type = typename Rotor::value_type;
		using index_type = typename Rotor::index_type;
		using size_type = typename Rotor::size_type;
		static auto immediately_add(Rotor&, const Rotor&) -> void;
		static auto immediately_multiply(Rotor&, const Rotor&) -> void;
		static auto immediately_divide(Rotor&, value_type) -> void;
		static auto minus(const Rotor&)->Rotor;
		static auto conjugate(const Rotor&)->Rotor;
		static auto calc_w(index_type, size_type)->Rotor;
	};

	template<class T>
	struct rotor_traits<std::complex<T>> {
		using value_type = T;
		using index_type = size_t;
		using size_type = size_t;

		static auto immediately_add(
			std::complex<T>& left, const std::complex<T>& right) -> void {
			left += right;
		}

		static auto immediately_multiply(
			std::complex<T>& left, const std::complex<T>& right) -> void {
			left *= right;
		}

		static auto immediately_divide(
			std::complex<T>& left, value_type right) -> void {
			left /= right;
		}

		static auto minus(const std::complex<T>& c) -> std::complex<T> {
			return -c;
		}

		static auto conjugate(const std::complex<T>& c) -> std::complex<T> {
			return std::conj(c);
		}

		static auto calc_w(index_type k, size_type n) -> std::complex<T> {
			return exp(std::complex<T>(0, (-2.0*3.14159*k) / static_cast<T>(n)));
		}
	};

	template<class SignalArray>
	struct signal_array_traits {
		using size_type = typename SignalArray::size_type;
		using value_type = typename SignalArray::value_type;
		using iterator_type = typename SignalArray::iterator;
		using index_type = size_t;

		static auto begin(SignalArray& sa)->iterator_type
		{
			return sa.begin();
		}

		static auto end(SignalArray& sa)->iterator_type
		{
			return sa.end();
		}

		static auto size(const SignalArray& sa)->size_type
		{
			return sa.size();
		}

		static SignalArray construct(size_type size)
		{
			return SignalArray(size);
		}

		static auto at(SignalArray& sa, index_type index)->value_type&
		{
			return sa.at(index);
		}

		static auto at(const SignalArray& sa, index_type index) -> const value_type&{
			return sa.at(index);
		}
	};

	template<class Rotor, class SignalArray>
	class fast_basic_transform {
	public:
		using r_traits = rotor_traits<Rotor>;
		using sa_traits = signal_array_traits<SignalArray>;

		fast_basic_transform(typename sa_traits::size_type sa_size) :
			n_(sa_size), n_level_(0), index_list_(1, 0) {
			typename sa_traits::size_type n = n_;
			while (n != 1) {
				++n_level_;
				n >>= 1;
			}
			init_memos();
		}

		auto transform(const SignalArray& sa) -> SignalArray {
			return transform_impl(sa, false);
		}

		auto inverse_transform(const SignalArray& sa) -> SignalArray {
			auto dst = transform_impl(sa, true);
			auto iter = sa_traits::begin(dst);
			while (iter != sa_traits::end(dst)) {
				r_traits::immediately_divide(*iter, n_);
				++iter;
			}
			return dst;
		}

	private:
		auto transform_impl(const SignalArray& sa, bool is_inverse) -> SignalArray {
			auto src = sa_traits::construct(sa_traits::size(sa));
			auto src_iter = sa_traits::begin(src);
			for (auto index : index_list_) {
				*src_iter = sa_traits::at(sa, index);
				++src_iter;
			}

			auto dst = sa_traits::construct(sa_traits::size(sa));
			for (typename sa_traits::index_type i = 0; i < n_level_; ++i) {
				dst = sa_traits::construct(n_);
				typename sa_traits::index_type block_num = (1 << (n_level_ - i - 1));
				for (typename sa_traits::index_type j = 0; j < block_num; ++j) {
					typename sa_traits::index_type wing = 1 << i;
					for (typename sa_traits::index_type k = 0; k < wing; ++k) {
						typename sa_traits::index_type index = k + 2 * wing*j;
						r_traits::immediately_add(
							sa_traits::at(dst, index), sa_traits::at(src, index));
						r_traits::immediately_add(
							sa_traits::at(dst, index + wing), sa_traits::at(src, index));

						auto w = is_inverse ?
							r_traits::conjugate(sa_traits::at(w_table_[i], k))
							: sa_traits::at(w_table_[i], k);
						r_traits::immediately_multiply(w, sa_traits::at(src, index + wing));
						r_traits::immediately_add(sa_traits::at(dst, index), w);
						r_traits::immediately_add(
							sa_traits::at(dst, index + wing), r_traits::minus(w));
					}
				}
				src = dst;
			}
			return dst;

		}

		auto init_memos() -> void {
			init_index_list();
			init_w_table();
		}

		auto init_index_list() -> void {
			index_list_.reserve(n_);
			for (unsigned int i = 0; i < n_level_; ++i) {
				index_list_.insert(index_list_.end(), index_list_.begin(), index_list_.end());
				for (unsigned int j = 0; j < index_list_.size(); ++j) {
					index_list_[j] *= 2;
				}
				for (unsigned int j = index_list_.size() / 2; j < index_list_.size(); ++j) {
					index_list_[j] += 1;
				}
			}
		}

		auto init_w_table() -> void {
			w_table_.clear();
			for (typename sa_traits::size_type i = 0; i < n_level_; ++i) {
				typename sa_traits::size_type wing = (1 << i);
				auto w_table_row = sa_traits::construct(wing);
				for (typename sa_traits::size_type k = 0; k < wing; ++k) {
					w_table_row[k] = r_traits::calc_w(k, 2 << i);
				}
				w_table_.push_back(w_table_row);
			}
		}

		typename sa_traits::size_type n_, n_level_;
		std::vector<typename sa_traits::index_type> index_list_;
		std::vector<SignalArray> w_table_;
	};
}
#endif

