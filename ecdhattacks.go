package dhpals

import (
	"bytes"
	"math/big"

	"github.com/dnkolegov/dhpals/elliptic"
	"github.com/dnkolegov/dhpals/x128"
)

var BigZero = big.NewInt(0)
var BigOne = big.NewInt(1)

func findSmallFactors(cofactor *big.Int) []*big.Int {
	var factors []*big.Int
	for i := 2; i < 1<<16; i++ {
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

			for a := BigOne; a.Cmp(order) <= 0; a.Add(a, BigOne) {

				cX, cY := curve.ScalarMult(x, y, a.Bytes())
				cur := append(cX.Bytes(), cY.Bytes()...)
				k := mixKey(cur)

				if bytes.Compare(msg, k) == 0 {
					A = append(A, a)
					N = append(N, order)
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

func runECDHTwistAttack(ecdh func(x *big.Int) []byte, getPublicKey func() (*big.Int, *big.Int), privateKeyOracle func(*big.Int) *big.Int) (priv *big.Int) {
	panic("not implemented")
	return
}

type twistPoint struct {
	order *big.Int
	point *big.Int
}

// findAllPointsOfPrimeOrderOnX128 finds a point with a specified order for u^3 + A*u^2 + u in GF(p).
func findAllPointsOfPrimeOrderOnX128() (points []twistPoint) {
	// It is known, that both curves contain 2*p+2 points: |E| + |T| = 2*p + 2
	panic("not implemented")
	x128.ScalarBaseMult(big.NewInt(1).Bytes())
	return
}

// catchKangarooOnMontgomeryCurve implements Pollard's kangaroo algorithm on a curve.
func catchKangarooOnCurve(curve elliptic.Curve, bx, by, x, y, a, b *big.Int) (m *big.Int, err error) {
	// k is calculated based on a formula in this paper: https://arxiv.org/pdf/0812.0789.pdf
	panic("not implemented")
	return
}
