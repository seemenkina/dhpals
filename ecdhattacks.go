package dhpals

import (
	"bytes"
	"context"
	"crypto/rand"
	"fmt"
	"github.com/davecgh/go-spew/spew"
	"math"
	"math/big"
	"runtime"
	"sync"

	"github.com/dnkolegov/dhpals/elliptic"
	"github.com/dnkolegov/dhpals/x128"
)

var (
	BigZero = big.NewInt(0)
	BigOne  = big.NewInt(1)
	BigTwo  = big.NewInt(2)
)

func findSmallFactors(cofactor *big.Int) []*big.Int {
	var factors []*big.Int

	for i := 3; i < 1<<16; i += 2 {
		f := new(big.Int).SetInt64(int64(i))
		tmp := new(big.Int).Mod(cofactor, f)

		if tmp.Cmp(BigZero) == 0 {
			factors = append(factors, f)
			for new(big.Int).Mod(cofactor, f).Cmp(BigZero) == 0 {
				cofactor = new(big.Int).Div(cofactor, f)
			}
		}
		if f.Cmp(cofactor) == 1 {
			break
		}
	}
	return factors
}

func findGenerator(order *big.Int, curve elliptic.Curve) (*big.Int, *big.Int) {

	for {
		x, y := elliptic.GeneratePoint(curve)

		gx, gy := curve.ScalarMult(x, y, new(big.Int).Div(curve.Params().N, order).Bytes())
		if gx.Cmp(BigZero) == 0 && gy.Cmp(BigZero) == 0 {
			continue
		} else {
			return gx, gy
		}
	}
}

func bruteHashICA(curve elliptic.Curve, x, y, posibleKey *big.Int) []byte {
	cX, cY := curve.ScalarMult(x, y, posibleKey.Bytes())
	cur := append(cX.Bytes(), cY.Bytes()...)
	currentMsg := mixKey(cur)
	return currentMsg
}

func runECDHInvalidCurveAttack(ecdh func(x, y *big.Int) []byte) (priv *big.Int) {

	badCurve := []elliptic.Curve{elliptic.P128V1(), elliptic.P128V2(), elliptic.P128V3()}
	var A, N []*big.Int

	for _, curve := range badCurve {
		factors := findSmallFactors(curve.Params().N)
		for _, order := range factors {

			x, y := findGenerator(order, curve)

			msg := ecdh(x, y)

			for possibleKey := BigOne; possibleKey.Cmp(order) <= 0; possibleKey = new(big.Int).Add(possibleKey, BigOne) {
				currentMsg := bruteHashICA(curve, x, y, possibleKey)

				if bytes.Equal(msg, currentMsg) {
					flag := true
					for i := 0; i < len(A) && i < len(N); i++ {
						if A[i].Cmp(possibleKey) == 0 && N[i].Cmp(order) == 0 {
							flag = false
							break
						}
					}
					if flag {
						A = append(A, possibleKey)
						N = append(N, order)
					}
				}

			}

		}

	}

	var tmpA []*big.Int
	var tmPN []*big.Int
	tmpA = append(tmpA, A[0])
	tmPN = append(tmPN, N[0])

	for i := 1; i < len(N); i++ {
		if N[i-1].Cmp(N[i]) == 0 {
			continue
		} else {
			tmpA = append(tmpA, A[i])
			tmPN = append(tmPN, N[i])
		}
	}

	x, _, err := crt(tmpA, tmPN)
	if err != nil {
		println(err)
		return nil
	}

	return x
}

func runECDHSmallSubgroupAttack(curve elliptic.Curve, ecdh func(x, y *big.Int) []byte) (priv *big.Int) {
	panic("not implemented")
	return
}

type twistPoint struct {
	order *big.Int
	point *big.Int
}

// findSmallTwistFactors finds all factors of r that are smaller than 1 << 24 exclusive and greater than 1 (if any)
func findSmallTwistFactors(r *big.Int) []*big.Int {
	var factors []*big.Int

	for i := 3; i < 1<<24; i += 2 { // finding all divisors of r
		f := new(big.Int).SetInt64(int64(i))
		tmp := new(big.Int).Mod(r, f)

		if tmp.Cmp(BigZero) == 0 {
			factors = append(factors, f)
			for new(big.Int).Mod(r, f).Cmp(BigZero) == 0 {
				r = new(big.Int).Div(r, f)
			}
		}
		if f.Cmp(r) == 1 {
			break
		}
	}

	return factors
}

// montgomeryPoly function calculates the value of point v = sqrt(u^3 + Au^2 + u) mod P
func montgomeryPoly(u *big.Int) *big.Int {
	// u^3
	u3 := new(big.Int).Mul(u, u)
	u3.Mul(u3, u)

	// A*u^2
	u2 := new(big.Int).Mul(u, u)
	u2.Mul(u2, new(big.Int).Set(x128.A))

	// (u^3 + Au^2 + u) mod P
	u3.Add(u3, u2)
	u3.Add(u3, u)
	u3.Mod(u3, x128.P)

	v := new(big.Int).ModSqrt(u3, new(big.Int).Set(x128.P))

	return v
}

// Finds a generator of a subgroup of E(GF(p)) of required order
// param order -- the order the generator must have, it must be a divisor of the order of the twist of this curve
// return the u coordinate of a generator of the twist curve satisfying the order given
func findTwistGenerator(order, twistOrder *big.Int) *big.Int {

	ord := new(big.Int).Div(twistOrder, order)

	for {
		u, err := rand.Int(rand.Reader, x128.P)
		if err != nil {
			continue
		}

		if v := montgomeryPoly(u); v == nil { // We are on the twist of this curve
			possibleGen := x128.ScalarMult(u, ord.Bytes())
			if possibleGen.Cmp(BigZero) == 0 {
				continue
			} else {
				return possibleGen
			}
		}
	}
}

// func getTwistOrder return order of the twist: 2 + 2*p - q
func getTwistOrder(p, q *big.Int) *big.Int {
	o := new(big.Int).Mul(p, BigTwo)
	o = o.Add(o, BigTwo)
	return o.Sub(o, q)
}

// findAllPointsOfPrimeOrderOnX128 finds a point with a specified order for u^3 + A*u^2 + u in GF(p).
func findAllPointsOfPrimeOrderOnX128() (points []twistPoint) {
	// It is known, that both curves contain 2*p+2 points: |E| + |T| = 2*p + 2
	twistOrder := getTwistOrder(x128.P, x128.N)
	factors := findSmallTwistFactors(twistOrder)
	spew.Dump(twistOrder, factors)
	for _, order := range factors {
		u := findTwistGenerator(order, twistOrder)
		points = append(points, twistPoint{order: order, point: u})
	}
	return
}

func getK(a, b *big.Int) *big.Int {
	t0 := new(big.Int).Sub(b, a)

	t0.Sqrt(t0)

	t1 := math.Log2(float64(t0.Int64()))
	t2 := math.Log2(t1)

	res := t1 + t2 - 2
	return new(big.Int).SetUint64(uint64(res))
}

func f(u, k *big.Int) *big.Int {
	// res = 2^(u % k)
	res := new(big.Int).Lsh(big.NewInt(1), uint(u.Uint64()%k.Uint64()))
	return res
}

func computeN(k *big.Int) *big.Int {
	N := big.NewInt(0)
	for i := big.NewInt(0); i.Cmp(k) <= 0; i = new(big.Int).Add(i, BigOne) {
		N.Add(N, f(i, k))
	}
	if new(big.Int).Rsh(k, 2).Cmp(Big0) != 0 {
		N.Div(N, new(big.Int).Rsh(k, 2))
	}

	return N
}

func MontgomeryToWeierstrass(u *big.Int) (*big.Int, *big.Int, error) {
	v := montgomeryPoly(u)
	if v == nil {
		return nil, nil, fmt.Errorf("no v coordinate on curve x128 for u coordinate")
	}
	//yS := []*big.Int{v, new(big.Int).Neg(v)}
	x := new(big.Int).Add(u, big.NewInt(178))
	if elliptic.P128().IsOnCurve(x, v) {
		fmt.Printf("V\n")
		return x, v, nil
	} else {
		if elliptic.P128().IsOnCurve(x, new(big.Int).Neg(v)) {
			fmt.Printf("-V\n")
			return x, new(big.Int).Neg(v), nil
		}
	}
	return nil, nil, fmt.Errorf("x, v are not on P128 curve")
}

// catchKangarooOnMontgomeryCurve implements Pollard's kangaroo algorithm on a curve.
func catchKangarooOnCurve(curve elliptic.Curve, bx, by, x, y, a, b *big.Int) (m *big.Int, err error) {
	// K is calculated based on a formula in this paper: https://arxiv.org/pdf/0812.0789.pdf
	K := getK(a, b)
	N := computeN(K)
	fmt.Printf("k: %d\nN: %d\n", K, N)
	// xTame, yTame := 0, b * base
	xTame := big.NewInt(0)
	xyTame, yyTame := curve.ScalarMult(bx, by, b.Bytes())

	for i := uint64(0); i < N.Uint64(); i++ {
		xTame.Add(xTame, f(xyTame, K))
		tmpX, tmpY := curve.ScalarMult(bx, by, f(xyTame, K).Bytes())

		xyTame, yyTame = curve.Add(xyTame, yyTame, tmpX, tmpY)
	}

	xB, yB := curve.ScalarMult(bx, by, new(big.Int).Add(b, xTame).Bytes())

	if xyTame.Cmp(xB) != 0 && yyTame.Cmp(yB) != 0 {
		return nil, fmt.Errorf("yTame == (b + xTame) * base should be true")
	}
	fmt.Printf("xTame: %d\nxyTame: %d\nyyTame: %d\n", xTame, xyTame, yyTame)
	// xWild, yWild := 0, (x,y)
	xWild := big.NewInt(0)
	xyWild, yyWild := x, y

	upperLimit := new(big.Int).Add(xTame, b)

	for xWild.Cmp(upperLimit) < 0 {

		xWild.Add(xWild, f(xyWild, K))

		tmpX, tmpY := curve.ScalarMult(bx, by, f(xyWild, K).Bytes())
		xyWild, yyWild = curve.Add(xyWild, yyWild, tmpX, tmpY)

		if xyTame.Cmp(xyWild) == 0 && yyTame.Cmp(yyWild) == 0 {
			fmt.Printf("xyWild: %d\nyyWild: %d\n", xyWild, yyWild)
			m = new(big.Int).Sub(xTame, xWild)
			m.Add(m, b)
			return m, nil
		}
	}

	return nil, fmt.Errorf("no result was found")
}

func runECDHTwistAttack(ecdh func(x *big.Int) []byte, getPublicKey func() *big.Int, privateKeyOracle func(*big.Int) *big.Int) (priv *big.Int) {

	mrw := sync.RWMutex{}
	var A, N []*big.Int
	points := findAllPointsOfPrimeOrderOnX128()

	type dataIn struct {
		cancel      context.CancelFunc
		possibleKey *big.Int
		msg         []byte
		point       twistPoint
	}

	nWorkers := runtime.NumCPU()
	chanIn := make(chan dataIn, nWorkers)

	var wg sync.WaitGroup
	wg.Add(nWorkers)

	for i := 0; i < nWorkers; i++ {
		go func() {
			defer wg.Done()

			for data := range chanIn {
				currentMsg := bruteHash(data.point, data.possibleKey)
				if bytes.Equal(data.msg, currentMsg) {
					mrw.RLock()
					flag := true
					for i := 0; i < len(A) && i < len(N); i++ {
						if A[i].Cmp(data.possibleKey) == 0 && N[i].Cmp(data.point.order) == 0 {
							flag = false
							break
						}
					}
					mrw.RUnlock()
					if flag {
						//fmt.Printf("point: %d, order %d\n", data.possibleKey, data.point.order)
						mrw.Lock()
						A = append(A, data.possibleKey)
						N = append(N, data.point.order)
						mrw.Unlock()
						data.cancel()
					}
				}

			}
		}()
	}

	for _, point := range points {
		msg := ecdh(point.point)
		ctx, cancel := context.WithCancel(context.Background())
		for posKey := BigOne; posKey.Cmp(point.order) <= 0; posKey = new(big.Int).Add(posKey, BigOne) {
			select {
			case <-ctx.Done():
				break
			case chanIn <- dataIn{
				possibleKey: posKey,
				msg:         msg,
				point:       point,
				cancel:      cancel,
			}:
				// pass
			}
		}
		cancel()
	}
	close(chanIn)
	wg.Wait()
	//for i := 0; i < len(A) && i < len(N); i++ {
	//	fmt.Printf("A: %d N:%d\n", A[i], N[i])
	//}

	var tmpA []*big.Int
	var tmPN []*big.Int
	tmpA = append(tmpA, A[0])
	tmPN = append(tmPN, N[0])

	for i := 1; i < len(N); i++ {
		if N[i-1].Cmp(N[i]) == 0 {
			continue
		} else {
			tmpA = append(tmpA, A[i])
			tmPN = append(tmPN, N[i])
		}
	}

	for i := 0; i < len(tmpA) && i < len(tmPN); i++ {
		fmt.Printf("A: %d N:%d\n", tmpA[i], tmPN[i])
	}
	candX, candN, err := crt(tmpA, tmPN)
	if err != nil {
		fmt.Print(fmt.Errorf("error: %v", err))
		return nil
	}
	fmt.Printf("x: %d\nn: %d\n", candX, candN)

	curve := elliptic.P128()

	baseNew := x128.ScalarBaseMult(candN.Bytes())
	baseNewx, baseNewy, err := MontgomeryToWeierstrass(baseNew)
	if err != nil {
		fmt.Print(fmt.Errorf("error: %v", err))
		return nil
	}
	fmt.Printf("xBase: %d\nyBase: %d\n", baseNewx, baseNewy)

	bobPubKey := getPublicKey()
	xNew, yNew, err := MontgomeryToWeierstrass(bobPubKey)
	if err != nil {
		fmt.Print(fmt.Errorf("error: %v", err))
		return nil
	}
	baseNx, baseNy := curve.ScalarBaseMult(new(big.Int).Neg(candX).Bytes())
	xNew, yNew = curve.Add(xNew, yNew, baseNx, baseNy)
	fmt.Printf("xNew: %d\nyNew: %d\n", xNew, yNew)

	a := big.NewInt(0)
	b := new(big.Int).Sub(x128.Q, big.NewInt(1))
	b.Div(b, candN)
	fmt.Printf("a: %d\nb: %d\n", a, b)

	m, err := catchKangarooOnCurve(curve, baseNx, baseNy, xNew, yNew, a, b)
	if err != nil {
		fmt.Print(fmt.Errorf("error: %v", err))
		return nil
	}

	priv = new(big.Int).Add(candX, new(big.Int).Mul(m, candN))
	return priv
}

func bruteHash(point twistPoint, possibleKey *big.Int) []byte {
	currentU := x128.ScalarMult(point.point, possibleKey.Bytes())
	currentMsg := mixKey(currentU.Bytes())
	return currentMsg
}
