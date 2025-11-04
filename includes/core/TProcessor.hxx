#pragma once

#include "libs.hh"
#include "AuxFunctions.hh"
#include <tuple>
#include <type_traits>
#include <utility>
#include "TContainer.hxx"

struct TProcessorBase {
	TProcessorBase() = default;
	TProcessorBase(const TProcessorBase& ) = default;
	TProcessorBase& operator=(const TProcessorBase& other) = default;
	TProcessorBase(TProcessorBase&& other) = default;
	TProcessorBase& operator=(TProcessorBase&& other) = default;

	/* virtual void ProcessEntry() noexcept = 0; */

	virtual ~TProcessorBase() = default;
};

template<typename T> struct TContainer;
template<typename T> struct TRawContainer;

template<class>
struct TProcessor; /* Undefined. */

/**
 * Base class of a single analysis subprocess.
 * Children will implement `ProcessEntry(...)` method to map the data
 * from combination of Input structures to the one unique Output structure.
 */
template<typename Out, typename... Ins>
struct TProcessor<Out(Ins...)> : TProcessorBase {
	static_assert((std::disjunction_v<
			util::is_base_of_template<TContainer, Ins>,
			util::is_base_of_template<TRawContainer, Ins>> && ...), 
		"Input type(s) must inherit from (or be) TContainer<T> / TRawContainer<T>.");
	static_assert(util::is_base_of_template<TContainer, Out>::value, "Output type must inherit from (or be) TContainer<T>.");

	Out out;
	std::tuple<Ins...> in;
	
	TProcessor() = default;

	/** 
	 * Only other acceptable param ctor is to take ownership from an existing output object.
	 * The output object doesn't need to be specifically `std::move`'ed to be less verbose in user code.
	 * But the Processor does take full ownership of the resource behind this reference.
	 * Input objects, might be shared between different subprocesses - here we simply make our own copy.
	 */
	explicit TProcessor(Out& _out, const Ins&... ins) : 
		out(std::move(_out)), in(std::make_tuple(ins...)) {}

	TProcessor(const TProcessor& rhs) : TProcessorBase(rhs), 
		out(rhs.out), 
		in(rhs.in) /* _vc from each input will just get copied - sharing the *same* pointers. */ 
		{
			out._vc.clear();

			/* https://en.cppreference.com/w/cpp/memory/shared_ptr/shared_ptr.html -- Case [13]
			 * template< class Y, class Deleter > shared_ptr( std::unique_ptr<Y, Deleter>&& r );*/
			for(const std::shared_ptr<TOnceBase>& v  : rhs.out._vc)
				out._vc.emplace_back( v->Clone() );

			/* Previous for-loop constructs the objects, this next dynamic dispatch
			 * will just give the raw pointer handles back to the user. */ 
			out.Setup();
		}
	/*  ^^^^^ Now, each `_vc` is completely unique in the output container. */
	
	/* Identical logic for copy-assignment op */
	TProcessor& operator=(const TProcessor& rhs) {
		this->in = rhs.in;
		this->out = rhs.out;
		this->out._vc.clear();
		
		for(const std::shared_ptr<TOnceBase>& v  : rhs.out._vc)
			out._vc.emplace_back( v->Clone() );
		
		out.Setup();
		return *this;
	}

	TProcessor(TProcessor&& )            noexcept = default;	
	TProcessor& operator=(TProcessor&& ) noexcept = default;	
	
	~TProcessor() = default;

	/**
	 * Returns a reference to the input container.
	 */
	template<u32 N = 0> 
	decltype(auto) GetInput() const& {
		static_assert(N < std::tuple_size_v<decltype(in)>, "Accessing N >= tuple size.");
		return ( std::get<N>(in).inner() );
	}

	template<u32 N = 0> 
	decltype(auto) GetInput() & {
		static_assert(N < std::tuple_size_v<decltype(in)>, "Accessing N >= tuple size.");
		return ( std::get<N>(in).inner() );
	}
	
	/* These calls are sent during the final collection. Stringent type checks
	 * are kept, as runtime isn't sacrificed too much. */
	void Collect(const TProcessor& rhs) {
		std::vector<std::shared_ptr<TOnceBase>>       & lvc = this->out._vc;
		const std::vector<std::shared_ptr<TOnceBase>> & rvc = rhs.out._vc;
		if( lvc.size() !=  rvc.size() )
			ERROR("(%s) trying to collect but output object named \'%s\' has unmatching sizes. %zu != %zu",
				_SELF_TYPE_CSTR, this->out.GetName(), lvc.size(), rvc.size());
		
		for(int i=0; i<(int)lvc.size(); ++i)
			lvc[i]->Collect( *rvc[i] );
	}

}; // TProcessor 
