/* Timer   - A simple timer class 
	 Created - 16/08/2016
*/

#ifndef TIMER_H
#define TIMER_H

#include <ctime>
#include <iostream>

namespace TSL
{
	// A simple timer class for timing methods
	
	class Timer
	{
		
	private:
		clock_t START_TIME;						// Start time
		clock_t PAUSED_TIME;					// Time that the clock was paused
		bool STARTED;									// Boolean stating if the clock has been started
		bool PAUSED;									// Boolean stating if the clock has been paused

	public:

		/* ----- Constructors and Destructor ----- */

		// Constructor
		Timer() : START_TIME( 0 ), PAUSED_TIME( 0 ), STARTED( false ), PAUSED( false )
		{} 

		// Destructor
		~Timer()
		{}

		/* ----- Methods ----- */

		bool is_started();						// Check if the timer has started
		bool is_stopped();						// Check if the timer has stopped
		bool is_paused();							// Check if the timer is paused	
		bool is_active();							// Check if the timer is still going

		void pause();									// Pause the timer
		void resume();								// Resume the timer
		void start();									// Start the timer
		void stop();									// Stop the timer
		void reset();									// Reset the timer
		void print() const;						// Output the time to the screen

		double get_time() const;			// Return the time in ms

		clock_t get_ticks();					// Return the number of clock ticks
	
	}; // End of Timer class

	/* ----- Inline definitions ----- */

	inline bool Timer::is_started()
	{
		return STARTED;
	}

	inline bool Timer::is_stopped()
	{
		return !STARTED;
	}

	inline bool Timer::is_paused()
	{
		return PAUSED;
	}

	inline bool Timer::is_active()
	{
		return !PAUSED & STARTED;
	}

	inline void Timer::pause()
	{
		if( PAUSED || !STARTED )
		{
			return;
		}
		PAUSED = true;
		PAUSED_TIME = clock();
	}

	inline void Timer::resume()
	{
		if( !PAUSED )
		{
			return;
		}
		PAUSED = false;
		START_TIME += clock() - PAUSED_TIME;
	}

	inline void Timer::start()
	{
		if( STARTED )
		{
			return;
		}
		STARTED = true;
		PAUSED = false;
		START_TIME = clock();
	}

	inline void Timer::stop()
	{
		STARTED = false;
	}

	inline void Timer::reset()
	{
		PAUSED = false;
		START_TIME = clock();
	}

	inline void Timer::print() const
	{
		std::cout.precision(4);
		Timer temp( *this );
		const double elapsed_time_in_ms( temp.get_time() );
		if ( elapsed_time_in_ms > 1000 )
        {
          std::cout << "  * TOTAL CPU time taken = " << elapsed_time_in_ms / 1000. << " s\n";
        }
        else
        {
         std::cout << "  * TOTAL CPU time taken = " << elapsed_time_in_ms << " ms\n";
        }	
	}
	
	inline double Timer::get_time() const
	{
		Timer temp( *this );
		return 1.e3 * temp.get_ticks() / CLOCKS_PER_SEC;
	}
	
	inline clock_t Timer::get_ticks()
	{
		if( !STARTED )
		{
			return 0;
		}

		if( PAUSED )
		{
			return PAUSED_TIME - START_TIME;
		}

		return clock() - START_TIME;
	}


}	// End of namespace TSL

#endif
