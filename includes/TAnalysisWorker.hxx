#pragma once

#include <atomic>
#include <thread>
#include <type_traits>
#include "Rtypes.h"
#include "libs.hh"

class TProcessor;

namespace util {
	template<typename T, typename = void>
	struct has_process_entry : std::false_type {};

	template<typename T>
	struct has_process_entry<T, std::void_t<decltype(std::declval<T&>().ProcessEntry())>> 
		: std::is_same<
			decltype(std::declval<T&>().ProcessEntry()),
			void
		> {};

	template<typename T>
	struct is_movable 
		: std::conjunction<
			std::is_move_assignable<T>,
			std::is_move_constructible<T>
		  > {};
}

template<typename T>
class TAnalysisWorker final : public T {
	static_assert(util::has_process_entry<T>::value, "Type <T> needs a `void ProcessEntry()` method implemented!");
	static_assert(std::is_move_constructible<T>::value,  "Type <T> needs a move ctor.");
	static_assert(std::is_move_assignable<T>::value, "Type <T> needs move assignment op.");
	static_assert(std::is_base_of<TProcessor, T>::value, "Type <T> must inherit from TProcessor!.");

	template<typename... Ts> friend class TAnalysisPool;

private:
	std::thread _thread{};
	std::atomic<bool> _stop{false};
	std::atomic<bool> _has_work{false};

public:
	using T::T;

	TAnalysisWorker() = default;
	TAnalysisWorker(const TAnalysisWorker& ) = delete;
	TAnalysisWorker& operator=(const TAnalysisWorker& ) = delete;

	TAnalysisWorker(TAnalysisWorker&& other) noexcept 
		: T(std::move(static_cast<T&>(other))),
		_thread(std::move(other._thread)),
		_stop(other._stop.load()),
		_has_work(other._has_work.load()) {}
	
	TAnalysisWorker& operator=(TAnalysisWorker&& other) noexcept {
		if(this != other) {
			T::operator=(std::move(static_cast<T&>(other)));

			_thread = std::move(other._thread);
			_stop.store(other._stop.load());
			_has_work.store(other._has_work.load());
		}
		return *this;
	}

	~TAnalysisWorker() { if(_thread.joinable()) _thread.join(); }

	bool IsDone() const noexcept { return ! _has_work.load(std::memory_order_acquire); }
		
	void Start() {
		if(_thread.joinable())
			ERROR("Trying to (re)start the worker but thread is in a joinable state?");
		if(_stop.load(std::memory_order_acquire)) 
			ERROR("Trying to start the worker but it's marked as stopped?");

		_thread = std::thread {
			[this]() {
				while(! _stop.load(std::memory_order_acquire)) {
					if(_has_work.load(std::memory_order_acquire)) {
						T::ProcessEntry();
						_has_work.store(false, std::memory_order_release);
					}
				}
			}
		};
	}

	/**
	 * Serialize single objects (non-TTree) parts of the worker to the currently open directory.
	 */
	Int_t Write() { return T::Write(); }
};
