function F = implied_debt(bg, bgp, bh, bhp, fut, tau)
global vars

	[curr,z] = outcomes_period(bg,bgp,bh,bhp,fut);
	% disp(curr(vars.T))
	% disp(curr(vars.N))
	tau_realized = curr(vars.T) / curr(vars.N);

	F = tau - tau_realized;
end
