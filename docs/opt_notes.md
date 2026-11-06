# cnamaste L-E-M optimization: objective, parameters and latents

Notes on `hmrf.run_core_inference` and the HMM fit it calls (2026-10-31).

**TL;DR:** cnamaste alternates two steps, and they do not score the same objective. Copy states are
averaged over (marginalized) by Baum-Welch when the emission parameters are fit on clone pseudobulks.
Clone labels are always a single best guess (hard MAP) from a Potts spatial model. Viterbi is not used
anywhere in this loop.

"L-E-M" is read here as L = clone labelling, E = forward-backward, M = BFGS parameter update. In the code,
each outer pass runs E and M together until the HMM fit converges, then runs L once:
[hmrf.py:527-650](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/hmrf.py#L527-L650).

## Notation

- spots `i`, bins `g`, clones `c`, copy states `k`
- per spot and bin: read count `y_ig`, baseline `e_ig` (`base_nb_mean`), B-allele count `b_ig` out of
  total allele count `n_ig`
- `l_i`: clone label of spot `i`
- `z_cg`: copy state of clone `c` at bin `g`
- `theta = {mu_k, p_k, alpha, tau}`: per-state read-depth rate and BAF, plus NB and beta-binomial
  dispersions

## 1. Overall objective

The model is a hidden Markov random field: one HMM along the genome per clone, plus a Potts prior over
neighbouring spots. Each pass has two blocks.

### E/M block (labels fixed)

Spots are summed into one pseudobulk per clone ([hmrf.py:593](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/hmrf.py#L593)).
The clones are stacked end to end along the genome, so they share `theta` but each keeps its own path
`z_c`. The fit maximizes the marginal likelihood:

```
L_EM(theta) = sum_c log sum_{z_c} p(z_c | T)
              * prod_g NB(Y_cg; E_cg * mu_{z_cg} / Z_c, alpha) * BB(B_cg; N_cg, p_{z_cg}, tau)
```

Stage 1 (BAF clones) runs with read depth zeroed and `params="sp"`; stage 2 (RDR refinement) runs with
`params="smp"` ([run_cnamaste.py:595](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/scripts/run_cnamaste.py#L595),
[:972](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/scripts/run_cnamaste.py#L972)).

### L block (`theta` and decoded paths `zhat` fixed)

Each spot is scored under each clone's decoded path
([hmrf.py:178](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/hmrf.py#L178)), and the labels maximize:

```
L_L(l)  = sum_i u_i(l_i) + beta * sum_{i~j} w_ij [l_i == l_j]
u_i(c)  = omega_i * sum_g log NB(y_ig; e_ig * mu_{zhat_cg} / Z_c, alpha)
          + sum_g log BB(b_ig; n_ig, p_{zhat_cg}, tau)
```

- `beta` is `hmrf.spatial_weight` (0.75 in `config.yaml`).
- `omega_i` is the spot's count of bins with allele reads divided by its count of bins with a baseline
  ([hmrf.py:158](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/hmrf.py#L158)).
- The energy is `potts.adjacency_potts_energy`
  ([potts.py](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/potts.py)).

### Is there one joint objective? NO

Four things break it:

1. **Pseudobulk vs per-spot scoring.** The E/M block scores clone sums; the L block scores each spot.
   - Under Poisson and binomial models the sum carries all the information, so both give the same
     likelihood for `mu` and `p`.
   - Under negative binomial and beta-binomial it does not. Summing per-spot NBs gives variance
     `sum m_i + alpha * sum m_i^2`, while the pseudobulk NB assumes `m + alpha * m^2` with
     `m = sum m_i`. So the fitted `alpha` and `tau` mean different things in the two blocks.
2. **Reweighting.** `omega_i` reweights the read-depth term in the L block only.
3. **Copy states.** The E/M block averages over `z`; the L block plugs in a fixed `zhat`.
4. **Changing clone count.** Merging and the clone-size floor (`min_spots_per_clone`) change `C` between
   passes.

So it is alternating optimization of two related objectives, with no guarantee that any one score goes
up every pass. That is why the loop stops when the labels stop changing (ARI against the previous labels
`>= hmrf.ari_tolerance`, which is 1.0, meaning identical labels) or at `max_iter_outer`. It then forces
one merge pass and one final fit ([hmrf.py:629](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/hmrf.py#L629),
[hmrf.py:652](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/hmrf.py#L652)).

> [!TIP]
> Compare with K-means (ISL Algorithm 12.2): assign spots to clones, then recompute each clone's
> "centroid", here an HMM copy-number profile. K-means always improves because both of its steps reduce
> the same sum of squares. This loop gives that up because its two blocks use different likelihoods.

## 2. Parameters vs latent variables

| Quantity | Role | How it is handled | Code |
|---|---|---|---|
| `mu_k`, `p_k` | parameter, one per state | point estimate by BFGS | [hmm_nophasing.py:908-938](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/hmm_nophasing.py#L908-L938) |
| `alpha`, `tau` | parameter | a single value each, shared by all states | [hmm_nophasing.py:695-698](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/hmm_nophasing.py#L695-L698) |
| transitions `T`, start probabilities | fixed | `t = 0.9999999`. `"s"` is in the parameter vector, but `cost_fn` ignores it, so the start probabilities stay at their initial values | [hmm_nophasing.py:909](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/hmm_nophasing.py#L909) |
| `Z_c` (per-clone library normalizer) | computed, not free | `log Z_c = logsumexp_g(log mu_{zhat_cg} + log lambda_g)`, using the previous pass's hard states | [hmm_nophasing.py:817](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/hmm_nophasing.py#L817), [:1057](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/hmm_nophasing.py#L1057) |
| `z_cg` copy state | latent | **averaged over** for fitting `theta`; **fixed to one value** for the L block: `zhat = argmax_k gamma_cg(k)` | [hmm.py:202-205](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/hmm.py#L202-L205), [hmrf.py:557](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/hmrf.py#L557) |
| `l_i` clone label | latent | **single best labelling only**: ICM under the repo's `config.yaml`, OpenGM alpha-expansion as the code default. Never averaged over (`posterior=None`) | [hmrf.py:319](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/hmrf.py#L319) |
| `beta` | hyperparameter | fixed | `config.yaml` |

### How E and M interleave

- **E:** forward-backward gives `gamma_cg(k) = p(z_cg = k | X_c, theta)`. No pairwise posteriors are
  computed, because `T` is fixed.
- **M:** BFGS minimizes `-sum gamma * log p(X | k, theta)`
  ([hmm_nophasing.py:938](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/hmm_nophasing.py#L938)), with the analytic gradient
  from `hmm_mstep_grad.em_gradient`.
- **Interleaving:** the E step runs inside the BFGS callback every second BFGS iteration
  ([hmm_nophasing.py:893](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/hmm_nophasing.py#L893)), so it is a generalized EM.
  - At convergence (`gamma` equals the posterior under `theta` and the gradient is zero), the point is a
    stationary point of `L_EM` (Fisher's identity).
  - Along the way, BFGS optimizes an objective that keeps changing, so the classic EM guarantee of steady
    improvement does not apply.

### Viterbi? NO

`zhat` is the most likely state at each bin separately (posterior-marginal decoding). Viterbi instead
finds the single most likely whole path. With `t` this close to 1 the two usually agree, but they are not
guaranteed to. Viterbi only appears later, in the integer allele-specific decoder
[hmm_ab_decoder.py](https://github.com/michaelJwilson/cna_maste/blob/42b96ff/python/cnamaste/hmm_ab_decoder.py).

### Overall

The copy states get a soft E step (Baum-Welch) for the parameters and are then fixed for labelling. The
labels are hard throughout, which makes the outer loop a classification-EM / ICM scheme. The soft
alternative, a mean-field average over labels, is not implemented.

## References

- MacKay, *Information Theory, Inference, and Learning Algorithms* (2003),
  https://www.inference.org.uk/itprnn/book.pdf
  - Ch. 20 and 22: hard vs soft K-means
  - Ch. 25: forward-backward, and decoding the most likely state per position vs the most likely whole
    path
  - Ch. 31: Ising/Potts models
  - Ch. 33: mean field as the soft alternative for labels
- James et al., *An Introduction to Statistical Learning* 2nd ed. (2021), Sec. 12.4.1, Algorithm 12.2
  (K-means), https://www.statlearning.com, doi:10.1007/978-1-0716-1418-1
- Dempster, Laird, Rubin (1977), EM, doi:10.1111/j.2517-6161.1977.tb01600.x
- Neal & Hinton (1998), generalized/incremental EM, doi:10.1007/978-94-011-5014-9_12
- Rabiner (1989), Baum-Welch and Viterbi, doi:10.1109/5.18626
- Besag (1986), ICM, doi:10.1111/j.2517-6161.1986.tb01412.x
- Boykov, Veksler, Zabih (2001), alpha-expansion, doi:10.1109/34.969114
- Celeux & Govaert (1992), classification EM, doi:10.1016/0167-9473(92)90042-E
- Ma et al. (2024), CalicoST, whose HMM/HMRF alternation this loop follows,
  doi:10.1038/s41592-024-02438-9
