/*************************************************************************
 *                                                                       *
 * polyjam, a polynomial solver generator for C++                        *
 * Copyright (C) 2015 Laurent Kneip, The Australian National University  *
 *                                                                       *
 * This program is free software: you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation, either version 3 of the License, or     *
 * (at your option) any later version.                                   *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. *
 *                                                                       *
 *************************************************************************/

#include <polyjam/math/GaussJordan.hpp>

#define PRECISION 0.0000000001

void
polyjam::math::customPrint( double value )
{
  std::cout << value << " ";
}

void
polyjam::math::customPrint( const core::Coefficient & value )
{
  value.print();
  std::cout << " ";
}

double
polyjam::math::customGetZero( double value )
{
  return 0.0;
}

double
polyjam::math::customGetOne( double value )
{
  return 1.0;
}

double
polyjam::math::customGetPrecision( double value )
{
  return PRECISION;
}

polyjam::core::Coefficient
polyjam::math::customGetZero( core::Coefficient value )
{
  return value.zero();
}

polyjam::core::Coefficient
polyjam::math::customGetOne( const core::Coefficient & value )
{
  return value.one();
}

polyjam::core::Coefficient
polyjam::math::customGetPrecision( const core::Coefficient & value )
{
  return core::Coefficient(PRECISION);
}
