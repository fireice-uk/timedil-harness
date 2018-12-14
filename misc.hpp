/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdint>
#include <vector>
#include <string>
#include <regex>
#include <array>
#include <list>

/** split a string in a vector of strings
 *
 * Based on Stack Overflow post:
 *   source: https://stackoverflow.com/a/28142357
 *   author: Marcin
 *   date: Jan 25 '15
 *
 * @param input string to split
 * @param regex separator between two elements
 */
std::vector< std::string > splitString(
    std::string const & input,
    std::string const & delimiter = ","
)
{
    std::regex re( delimiter );
    // passing -1 as the submatch index parameter performs splitting
    std::sregex_token_iterator first{
        input.begin(),
        input.end(),
        re,
        -1
    };
    std::sregex_token_iterator last;

    return {
        first,
        last
    };
}

 struct Interval
{
    /** time slice configuration
     *
     * 0 = begin of the interval
     * 1 = end of the interval
     * 2 = period
     */
    std::array< uint32_t, 3 > values;

    std::string toString() const
    {
        std::string result;
        result = std::to_string(values[0]) + ":" +
            std::to_string(values[1]) + ":" +
            std::to_string(values[2]);
        return result;
    }

    /** set the value
     *
     * if str is empty the default value for the given index is selected
     *
     * @param idx index to set, range [0,3)
     * @param str value to set, can be empty
     */
    void setValue(uint32_t const idx, std::string const & str)
    {
        if(!str.empty())
        {
            uint32_t value = std::stoul( str );
            values.at( idx )  = value;
        }
    }

    //! create a interval slice instance
    Interval() :
        /* default: start:end:period
         * -1 stored as unsigned is the highest available unsigned integer
         */
        values( { 0, uint32_t( -1 ), 1 } )
    { }
};

/** check if a given iteration is in the interval list
 *
 * @param seqIntervals vector with intervals
 * @param iteration step to check
 * @return true if step is included in the interval list else false
 */
bool containsStep(
    std::vector< Interval > const & seqIntervals,
    uint32_t const iteration
)
{
    for(auto const & slice : seqIntervals)
    {
        if(
            iteration >= slice.values[ 0 ] &&
            iteration <= slice.values[ 1 ]
        )
        {
            uint32_t const timeRelativeToStart = iteration - slice.values[ 0 ];
            if( timeRelativeToStart % slice.values[ 2 ] == 0 )
                return true;
        }
    }
    return false;
}

/** check if string contains only digits
 *
 * @param str string to check
 * @return true if str contains only digits else false
 */
bool is_number( std::string const & str )
{
    return std::all_of(
        str.begin(),
        str.end(),
        ::isdigit
    );
}

/** removes all spaces in a string */
std::string removeSpaces( std::string value )
{
    value.erase(
	std::remove(
	    value.begin(),
	    value.end(),
	    ' '
	),
	value.end()
    );

    return value;
}


/** create a Interval out of an string
 *
 * Parse a comma separated list of time slices and creates a vector of Intervals.
 * Interval Syntax:
 *   - `start:stop:period`
 *   - a number ``N is equal to `::N`
 *
 * - start:stop:period means the interval [start,stop] and each period's iteration
 */
std::vector< Interval > toInterval( std::string const & str )
{
    std::vector< Interval > result;
    auto const seqOfSlices = splitString(
        str,
        ","
    );
    for( auto const & slice : seqOfSlices )
    {
        auto const sliceComponents = splitString(
            slice,
            ":"
        );


        // id of the component
        size_t n = 0;
        bool const hasOnlyPeriod = sliceComponents.size() == 1u;

		Interval interval;
        for( auto& component : sliceComponents )
        {
			interval.setValue(
                hasOnlyPeriod ? 2 : n,
                component
            );
            n++;
        }
        result.push_back( interval );

    }
    return result;
}


struct HashPowerFunction
{

    virtual void run(uint64_t & run, const std::string & name, const size_t step) = 0;

    virtual std::string getName() const = 0;

};

struct TimestampFunction
{

    virtual void run(uint64_t & fakeTime, const uint64_t solveTime, const std::string & name, const Interval& slice, const size_t step) = 0;

    virtual std::string getName() const = 0;

};


class FunctorConnector
{
private:
    using SeqOfIntervals = std::vector< Interval >;
    using PluginPair = std::pair<
        std::string,
        SeqOfIntervals
    >;
    using NotificationList = std::list< PluginPair >;

public:

    /**
     * Notifies plugins that data should be dumped.
     *
     * @param currentStep current simulation iteration step
     */
    void hashPower(const std::vector<HashPowerFunction*>& functors, uint32_t currentStep, uint64_t& hr)
    {

        for (NotificationList::iterator iter = notificationList.begin();
                iter != notificationList.end(); ++iter)
        {
            for( auto & f : functors)
            {
                const std::string fName = f->getName();
                const auto funcValue = splitString(iter->first, ",");
                if(
                    fName == funcValue[0] &&
                    containsStep(
                        iter->second,
                        currentStep
                    )
                )
                {
                    f->run(hr, iter->first, currentStep);
                }
            }
        }
    }

    void timeStamp(const std::vector<TimestampFunction*>& functors, uint32_t currentStep, uint64_t & fakeTime, const size_t solveTime)
    {

        for (NotificationList::iterator iter = notificationList.begin();
                iter != notificationList.end(); ++iter)
        {
            for( auto & f : functors)
            {
                const std::string fName = f->getName();
                const auto funcValue = splitString(iter->first, ",");
                if(
                    fName == funcValue[0] &&
                    containsStep(
                        iter->second,
                        currentStep
                    )
                )
                {
                    f->run(fakeTime, solveTime, iter->first, iter->second[0], currentStep);
                }
            }
        }
    }

    /** Set the notification period
     *
     * @param notifiedObj the object to notify, e.g. an IPlugin instance
     * @param period notification period
     */
    void addRules(std::string const & all_rules)
    {
	    const auto rules_no_spaces = removeSpaces(all_rules);
	    if(rules_no_spaces.empty())
		    return;

        auto single_rules = splitString(rules_no_spaces, ";");
        if( !single_rules.empty() )
        {

            for(const auto & r : single_rules)
            {
                auto operation = splitString(r, "\\|");

                SeqOfIntervals seqIntervals = toInterval( operation[1] );
                notificationList.push_back( std::make_pair(
                    operation[0],
                    seqIntervals
                ) );
            }
        }

    }


    static FunctorConnector& getInstance()
    {
        static FunctorConnector instance;
        return instance;
    }

	private:

    FunctorConnector()
    {

    }

    virtual ~FunctorConnector()
    {

    }

    NotificationList notificationList;
};
