package dhpals

import (
	"bytes"
	"crypto/rand"
	"fmt"
	"math/big"

	"github.com/davecgh/go-spew/spew"
	"github.com/dnkolegov/dhpals/elliptic"
	"github.com/dnkolegov/dhpals/x128"
)

var BigZero = big.NewInt(0)
var BigOne = big.NewInt(1)

func findSmallFactors(cofactor *big.Int) []*big.Int {
	var factors []*big.Int
	for i := 3; i < 1<<16; i += 2 {
		I := new(big.Int).SetInt64(int64(i))
		z := new(big.Int).Mod(cofactor, I)
		if z.Cmp(BigZero) == 0 {
			factors = append(factors, I)
			for new(big.Int).Mod(cofactor, I).Cmp(BigZero) == 0 {
				cofactor = new(big.Int).Div(cofactor, I)
			}
		}
		if I.Cmp(cofactor) == 1 {
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

func runECDHInvalidCurveAttack(ecdh func(x, y *big.Int) []byte) (priv *big.Int) {

	badCurve := []elliptic.Curve{elliptic.P128V1(), elliptic.P128V2(), elliptic.P128V3()}
	var A, N []*big.Int

	for _, curve := range badCurve {
		factors := findSmallFactors(curve.Params().N)
		for _, order := range factors {

			x, y := findGenerator(order, curve)

			msg := ecdh(x, y)

			for a := BigOne; a.Cmp(order) <= 0; a = new(big.Int).Add(a, BigOne) {
				cX, cY := curve.ScalarMult(x, y, a.Bytes())
				cur := append(cX.Bytes(), cY.Bytes()...)
				k := mixKey(cur)

				if bytes.Equal(msg, k) {
					flag := true
					for i := 0; i < len(A) && i < len(N); i++ {
						if A[i].Cmp(a) == 0 && N[i].Cmp(order) == 0 {
							flag = false
						}
					}
					if flag {
						A = append(A, a)
						N = append(N, order)
						break
					}
				}

			}

		}

	}

	x, _, err := crt(A, N)
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

func findSmallTwistFactors(factor *big.Int) []*big.Int {
	var factors []*big.Int
	for i := 3; i < 1<<24; i += 2 {
		I := new(big.Int).SetInt64(int64(i))
		z := new(big.Int).Mod(factor, I)
		if z.Cmp(BigZero) == 0 {
			factors = append(factors, I)
			for new(big.Int).Mod(factor, I).Cmp(BigZero) == 0 {
				factor = new(big.Int).Div(factor, I)
			}
		}
		if I.Cmp(factor) == 1 {
			break
		}
	}
	return factors
}

func findTwistGenerator(order *big.Int) *big.Int {

	twistOrder := getTwistOrder(x128.P, x128.N)
	ord := new(big.Int).Div(twistOrder, order)

	for {
		u, err := rand.Int(rand.Reader, x128.P)
		if err != nil {
			continue
		}
		u3 := new(big.Int).Mul(u, u)
		u3.Mul(u3, u)

		u2 := new(big.Int).Mul(u, u)
		u2.Mul(u2, x128.A)

		u3.Add(u3, u2)
		u3.Add(u3, u)
		u3.Mod(u3, x128.P)

		v := new(big.Int).ModSqrt(u3, x128.P)
		if v == nil {
			possibleGen := x128.ScalarMult(u, ord.Bytes())
			if possibleGen.Cmp(BigZero) == 0 {
				continue
			} else {
				return possibleGen
			}
		}
	}
}

func getTwistOrder(p, q *big.Int) *big.Int {
	o := new(big.Int).Add(big.NewInt(2), new(big.Int).Mul(big.NewInt(2), p))
	return o.Sub(o, q)
}

// findAllPointsOfPrimeOrderOnX128 finds a point with a specified order for u^3 + A*u^2 + u in GF(p).
func findAllPointsOfPrimeOrderOnX128() (points []twistPoint) {
	// It is known, that both curves contain 2*p+2 points: |E| + |T| = 2*p + 2
	factor := getTwistOrder(x128.P, x128.N)
	factors := findSmallTwistFactors(factor)

	for _, order := range factors {

		x := findTwistGenerator(order)
		points = append(points, twistPoint{order: order, point: x})
	}
	return
}

// catchKangarooOnMontgomeryCurve implements Pollard's kangaroo algorithm on a curve.
func catchKangarooOnCurve(curve elliptic.Curve, bx, by, x, y, a, b *big.Int) (m *big.Int, err error) {
	// k is calculated based on a formula in this paper: https://arxiv.org/pdf/0812.0789.pdf
	panic("not implemented")
	return
}

func runECDHTwistAttack(ecdh func(x *big.Int) []byte, getPublicKey func() (*big.Int, *big.Int), privateKeyOracle func(*big.Int) *big.Int) (priv *big.Int) {

	var A, N []*big.Int
	points := findAllPointsOfPrimeOrderOnX128()

	for _, point := range points {
		msg := ecdh(point.point)

		for a := BigOne; a.Cmp(point.order) <= 0; a = new(big.Int).Add(a, BigOne) {
			cU := x128.ScalarMult(point.point, a.Bytes())
			cMsg := mixKey(cU.Bytes())

			if bytes.Equal(msg, cMsg) {
				flag := true
				for i := 0; i < len(A) && i < len(N); i++ {
					if A[i].Cmp(a) == 0 && N[i].Cmp(point.order) == 0 {
						flag = false
					}
				}
				if flag {
					A = append(A, a)
					N = append(N, point.order)
					break
				}
			}

		}
	}

	spew.Dump(A, N)

	x, _, err := crt(A, N)
	if err != nil {
		println(err)
	}
	fmt.Printf("x:%d\n", x)

	return
}
