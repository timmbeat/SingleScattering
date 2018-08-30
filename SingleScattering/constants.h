#pragma once


namespace mcss
{
	
	using Real = long double;


	template<typename Real> constexpr Real pi()
	{
		return static_cast<Real>(3.14159265358979323846);
	}
	template<typename Real> constexpr Real k0low(Real alpha)
	{
		return static_cast<Real> (1 - 2 * exp(-2 / alpha)*(1 + ((4 - alpha) / alpha)*exp(-2 / alpha) + 
			(((24 - 12 * alpha + alpha * alpha) / (alpha*alpha))*exp(-4 / alpha))) + 
								 (((512 - 384 * alpha + 72 * alpha*alpha - 3 * alpha*alpha*alpha) / (alpha*alpha*alpha))*exp(-6 / alpha)));
	}

	template<typename Real> constexpr Real k0high(Real alpha)
	{
		Real const omalpha = (1 - alpha);
		Real const qomalpha = (1 - alpha)*(1 - alpha);
		return static_cast<Real> (sqrt(3 * omalpha)*(1-((2/5)*omalpha) - ((12/175)*qomalpha) - ((2/125)*qomalpha*omalpha)+((166/67375)*qomalpha*qomalpha)));
	}

	template<typename Real> constexpr Real v0(Real alpha)
	{
		auto const alph = alpha > 0.9999 ? 0.9999 : alpha;
		return alph > 0.56 ? static_cast<Real>(1/k0high<Real>(alph)) : static_cast<Real>(1/k0low<Real>(alph));
	}


}
