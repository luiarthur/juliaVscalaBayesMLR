// usage: sbt run n k (with intercept)
object mlr {
  def main(args: Array[String]) = {
  // set variables
    val (n,k) = (args(0).toInt, args(1).toInt) 

    def timer[R](block: => R): R = {  
        val t0 = System.nanoTime()
        val result = block    // call-by-name
        val t1 = System.nanoTime()
        println(Console.GREEN + "\nElapsed time: " + (t1 - t0) / 1E9 + "s\n" +
          Console.RESET)
        result
    }

    import breeze.linalg._
    import breeze.linalg.{DenseMatrix=>dmat,DenseVector=>dvec}
    import breeze.numerics._
    import breeze.stats.distributions._ // 1 draw: Dist.draw. n draws: Dist.sample(n)
    import breeze.stats.mean
    val Rand = scala.util.Random

    class State(val beta: dvec[Double], val s2: Double)
    val G = Gaussian(0,1)
    val U = Uniform(0,1)

    // generate data:
    val X = dmat.horzcat(
      dmat.ones[Double](n,1),
      dmat.fill(n,k-1)(Rand.nextGaussian)
    )
    val betaTrue = dvec.tabulate(k)(j => j + 1.0)
    val y = X * betaTrue + dvec.fill(n)(Rand.nextGaussian)

    val Xt = X.t
    val XXi = inv(Xt*X)
    val s2 = 10.0
    val a = 1.0
    val b = 1.0
    val B = 100000
    val keep = B / 10
    val burn = B - keep
    val csb = XXi(::,*) * 4.0
    val S = cholesky(csb)
    val css = 1.0
    var accb, accs = 0

    def ll(be: dvec[Double], sig2: Double): Double = {
      val c = y - X * be
      (c.t*c/sig2 + n*log(sig2)) / -2.0
    }

    def lpb(be: dvec[Double]): Double = be.t*XXi*be / (-2.0*s2)
    def lps(sig2: Double): Double = (a-1)*log(sig2) - sig2/b
    def mvrnorm(m: dvec[Double]) = m + S*G.samplesVector(k)

    // Update Functions:
    def update_State(s: State): State = {

      // Update beta
      val b_cand = mvrnorm(s.beta)
      val q = ll(b_cand,s.s2)+lpb(b_cand) - ll(s.beta,s.s2)-lpb(s.beta)
      val beta_new = if (q > log(U.draw)) { accb += 1; b_cand} else s.beta

      // Update s2
      val s2_cand= G.draw * sqrt(css) + s.s2
      val s2_new = if (s2_cand> 0) {
        val q = ll(beta_new, s2_cand)+lps(s2_cand) - ll(beta_new, s.s2)-lps(s.s2)
        if (q > log(U.draw)) {accs += 1; s2_cand} else s.s2
      } else s.s2

      new State( beta_new, s2_new )
    }

    def mh(i: Int, S: List[State]): List[State] = {
      if (i < B) {
        if (i % (B/10) == 0) print("\rProgress: "+round(i*100.0/B)+"%")
        val new_State = update_State(S.head) :: S
        mh(i+1, new_State)
      } else S.dropRight(burn)
    }

    val out = timer { mh(1, List(new State(dvec.zeros[Double](k), 1.0))) }
    val (beta_post, sig2_post) = out.map(s => (s.beta, s.s2)).unzip

    println("Acceptance beta: " + 100.0 * accb/B+"%")
    println("Acceptance sig2: " + 100.0 * accs/B+"%\n")
    println("Posterior mean sig2: " + sig2_post.map(x => x / keep.toDouble).reduce(_+_))
    println("Posterior mean beta:" )
    (beta_post.map(x => x / keep.toDouble).reduce(_+_)).foreach(println)
  }
}
