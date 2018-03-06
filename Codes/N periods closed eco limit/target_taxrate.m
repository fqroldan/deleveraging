function [F,z] = target_taxrate(tau,bgp,bh,bhp,fut)

	% disp(implied_debt(4,bgp,bh,bhp,fut,tau))
	bg = fzero( @(Bg) implied_debt(Bg,bgp,bh,bhp,fut,tau), bgp );

	[F,z] = outcomes_period(bg,bgp,bh,bhp,fut);
end
