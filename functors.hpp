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


namespace functors
{
	/** set the hash rate to X
	 *
	 * usage: "hrSet,50000"
	 */
	struct hrSet : public HashPowerFunction
	{
	    virtual void run(uint64_t & hr, const std::string & name, const size_t step)
	    {
		auto funcValue = splitString(name, ",");
		hr = std::stoull(funcValue[1]);
		std::cerr<<step<<": hrSet "<<hr<<std::endl;
	    }

	    virtual std::string getName() const
	    {
		return "hrSet";
	    }
	};

	/** multiply current hash rate with X
	 *
	 * usage: "hrMul,0.5"
	 *        "hrMul,1.5"
	 */
	struct hrMul : public HashPowerFunction
	{
	    virtual void run(uint64_t & hr, const std::string & name, const size_t step)
	    {
		auto funcValue = splitString(name, ",");
		hr *= std::stod(funcValue[1]);
		std::cerr<<step<<": hrMul "<<hr<<std::endl;
	    }

	    virtual std::string getName() const
	    {
		return "hrMul";
	    }
	};

	/** add an offset multiplied by a scaling factor to the block timestamp
	 *
	 * scaling_factor is `current_step - start_of_functor_interval`
	 *
	 * calculate: current_timestamp += scaling_factor * X
	 *
	 * usage:  addScaled,-1
	 */
	struct addScaled : public TimestampFunction
	{

	    virtual void run(uint64_t & fakeTime, const uint64_t solveTime, const std::string & name, const Interval& slice, const size_t step)
	    {
		auto funcValue = splitString(name, ",");
		int64_t offset = (step - slice.values[0] + 1) * std::stoll(funcValue[1]);
		std::cerr<<step<<": addScaled offset "<<offset<<std::endl;
		fakeTime += offset;
	    }

	    virtual std::string getName() const
	    {
		return "addScaled";
	    }

	};

	/** add an offset to the block timestamp
	 *
	 * calculate: current_timestamp += X
	 *
	 * usage:  add,100
	 */
	struct add : public TimestampFunction
	{

	    virtual void run(uint64_t & fakeTime, const uint64_t solveTime, const std::string & name, const Interval& slice, const size_t step)
	    {
		auto funcValue = splitString(name, ",");
		int64_t offset = std::stoll(funcValue[1]);
		std::cerr<<step<<": add offset "<<offset<<std::endl;
		fakeTime += offset;
	    }

	    virtual std::string getName() const
	    {
		return "add";
	    }

	};

	/** set a fixed block timestamp
	 *
	 * This functor allows to store the current block time and set it in the next iterations to the stored timestamp.
	 * This functor is a singleton and each time it is used in a independent rule it will use the timestamp from a previous store command.
	 *
	 * calculate: current_time = stored_timestamp
	 *
	 * usage:  set,s  -> stores the current block timestamp
	 *         set    -> without argument the last stored timestamp will be stet
	 */
	struct set : public TimestampFunction
	{

	    virtual void run(uint64_t & fakeTime, const uint64_t solveTime, const std::string & name, const Interval& slice, const size_t step)
	    {
		auto funcValue = splitString(name, ",");
		if(funcValue.size() == 2u && funcValue[1] == "s")
		{
		    std::cerr<<step<<": set store "<<fakeTime<<std::endl;
		    value = fakeTime;
		}
		else
		{
		    std::cerr<<step<<": set time to "<<value<<std::endl;
		    fakeTime = value;
		}
	    }

	    virtual std::string getName() const
	    {
		return "set";
	    }

	    uint64_t value = 0;
	};

}