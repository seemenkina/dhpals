// Package x128 implements the insecure Montgomery curve x128 defined in the Cryptopals challange 59.
package x128

import (
	"crypto/rand"
	"io"
	"math/big"
)

var (
	// A  - the a parameter.
	A = big.NewInt(534)
	// N - the order of the base point.
	N, _ = new(big.Int).SetString("233970423115425145498902418297807005944", 10)
	// P - the order of the underlying field.
	P, _ = new(big.Int).SetString("233970423115425145524320034830162017933", 10)
	// Q - the order of the subgroup.
	Q, _ = new(big.Int).SetString("29246302889428143187362802287225875743", 10)
	// U - the base point coordinate.
	U = big.NewInt(4)
	// V - the base point coordinate.
	V, _ = new(big.Int).SetString("85518893674295321206118380980485522083", 10)
	zero = big.NewInt(0)
	one  = big.NewInt(1)
	two  = big.NewInt(2)
	four = big.NewInt(4)
)

func ScalarBaseMult(k []byte) *big.Int {
	return ScalarMult(U, k)
}

func ScalarMult(in *big.Int, k []byte) *big.Int {
	return ladder(in, new(big.Int).SetBytes(k))
}

// polynomial returns u^3 + Au^2 + u
func polynomial(u *big.Int) *big.Int {
	u3 := new(big.Int).Mul(u, u)
	u3.Mul(u3, u)

	u2 := new(big.Int).Mul(u, u)
	u2.Mul(u2, A)

	u3.Add(u3, u2)
	u3.Add(u3, u)
	u3.Mod(u3, P)

	return u3
}

func IsOnCurve(u, v *big.Int) bool {
	// v^2 = u^3 + Au^2 + u
	newV := new(big.Int).Mul(v, v)
	newV.Mod(newV, P)

	return newV.Cmp(polynomial(u)) == 0
}

func cswap(x, y *big.Int, b bool) (u, v *big.Int) {
	if b {
		x, y = y, x
	}
	return x, y
}

func ladder(u, k *big.Int) *big.Int {
	u2, w2 := big.NewInt(1), big.NewInt(0)
	u3, w3 := new(big.Int).Set(u), big.NewInt(1)

	for i := P.BitLen() - 1; i >= 0; i-- {

		// b := 1 & (k >> i)
		b := new(big.Int).And(big.NewInt(1), new(big.Int).Rsh(k, uint(i)))

		// u2, u3 := cswap(u2, u3, b)
		// w2, w3 := cswap(w2, w3, b)
		u2, u3 = cswap(u2, u3, b.Cmp(big.NewInt(1)) == 0)
		w2, w3 = cswap(w2, w3, b.Cmp(big.NewInt(1)) == 0)

		// u3 := (u2*u3 - w2*w3)^2
		// w3 := u * (u2*w3 - w2*u3)^2

		// t := u2*u3 - w2*w3
		// tt := u2*w3 - w2*u3
		t := new(big.Int).Sub(new(big.Int).Mul(u2, u3), new(big.Int).Mul(w2, w3))
		tt := new(big.Int).Sub(new(big.Int).Mul(u2, w3), new(big.Int).Mul(w2, u3))

		u3 = new(big.Int).Mul(t, t)
		u3.Mod(u3, P)

		w3 = new(big.Int).Mul(tt, tt)
		w3 = new(big.Int).Mul(w3, u)
		w3.Mod(w3, P)

		// u2 := (u2^2 - w2^2)^2
		// w2 := 4*u2*w2 * (u2^2 + A*u2*w2 + w2^2)

		// t := u2^2
		// tt := w2^2
		// ttt := u2*w2

		t = new(big.Int).Mul(u2, u2)
		tt = new(big.Int).Mul(w2, w2)
		ttt := new(big.Int).Mul(u2, w2)

		u2.Sub(t, tt)
		u2.Mul(u2, u2)
		u2.Mod(u2, P)

		w2.Mul(ttt, four)

		// tttt := u2^2 + A*u2*w2 + w2^2 = t + A*ttt + tt
		tttt := new(big.Int).Mul(ttt, A)
		tttt.Add(tttt, tt)
		tttt.Add(tttt, t)

		w2.Mul(w2, tttt)
		w2.Mod(w2, P)

		// u2, u3 := cswap(u2, u3, b)
		// w2, w3 := cswap(w2, w3, b)
		u2, u3 = cswap(u2, u3, b.Cmp(big.NewInt(1)) == 0)
		w2, w3 = cswap(w2, w3, b.Cmp(big.NewInt(1)) == 0)
	}

	// return u2 * w2^(p-2)
	exp := new(big.Int).Exp(w2, new(big.Int).Sub(P, two), P)
	u2.Mul(u2, exp)
	u2.Mod(u2, P)
	return u2
}

func GenerateKey(rng io.Reader) (priv []byte, pub *big.Int, err error) {
	if rng == nil {
		rng = rand.Reader
	}

	bitSize := Q.BitLen()
	byteLen := (bitSize + 7) >> 3
	priv = make([]byte, byteLen)

	for pub == nil {
		_, err = io.ReadFull(rng, priv)
		if err != nil {
			return
		}
		if new(big.Int).SetBytes(priv).Cmp(Q) >= 0 {
			continue
		}

		pub = ScalarBaseMult(priv)
	}
	return
}
