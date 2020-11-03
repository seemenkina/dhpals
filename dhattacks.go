package dhpals

import (
	"bytes"
	"crypto/rand"
	"fmt"
	"math"
	"math/big"

	"github.com/davecgh/go-spew/spew"
)

func findSmallFactors(cofactor *big.Int) []*big.Int {
	var factors []*big.Int
	for i := 2; i < 1<<16; i++ {
		I := new(big.Int).SetInt64(int64(i))
		z := new(big.Int).Mod(cofactor, I)
		if z.Cmp(Big0) == 0 {
			factors = append(factors, I)
			for new(big.Int).Mod(cofactor, I).Cmp(Big0) == 0 {
				cofactor = new(big.Int).Div(cofactor, I)
			}
		}
		if I.Cmp(cofactor) == 1 {
			break
		}
	}
	return factors
}

// find an element of order order mod p
func findElementOfOrder(order, p *big.Int) (*big.Int, error) {
	h := Big1

	for h.Cmp(Big1) == 0 {
		buf, err := rand.Int(rand.Reader, p)
		if err != nil {
			return nil, err
		}
		if buf.Cmp(Big0) == 0 {
			buf = buf.Add(buf, Big1)
		}
		expBuf := new(big.Int).Sub(p, Big1)
		exp := new(big.Int).Div(expBuf, order)

		h = new(big.Int).Exp(buf, exp, p)
	}

	return h, nil
}

func bruteForcedMAC(digest []byte, h, order, p *big.Int) (*big.Int, error) {

	if new(big.Int).Exp(h, order, p).Cmp(Big1) != 0 {
		return nil, fmt.Errorf("element h should have had order order")
	}

	for i := int64(1); i <= order.Int64(); i++ {
		public := new(big.Int).Exp(h, new(big.Int).SetInt64(i), p)
		candidate := mixKey(public.Bytes())

		if bytes.Compare(candidate, digest) == 0 {
			return new(big.Int).SetInt64(i), nil
		}
	}

	return nil, fmt.Errorf(" could not find h^x such that h^x = K mod p ")
}

func runDHSmallSubgroupAttack(p, cofactor *big.Int, dh func(*big.Int) []byte) (priv *big.Int) {
	var B, R []*big.Int

	factors := findSmallFactors(cofactor)

	for _, order := range factors {
		h, err := findElementOfOrder(order, p)
		if err != nil {
			break
		}
		digest := dh(h)

		b, err := bruteForcedMAC(digest, h, order, p)
		if err != nil {
			fmt.Printf("Error: %v", err)
			break
		}
		B = append(B, b)
		R = append(R, order)
	}

	priv, _, err := crt(B, R)
	if err != nil {
		fmt.Printf("Error: %v", err)
		return nil
	}

	return priv
}

func computeN(k *big.Int) *big.Int {

	N := Big0
	for i := uint64(0); i < k.Uint64(); i++ {
		N = new(big.Int).Add(N, f(big.NewInt(int64(i)), k))
	}
	N = new(big.Int).Mul(big.NewInt(4), new(big.Int).Div(N, k))

	return N
}

func f(y, k *big.Int) *big.Int {
	return new(big.Int).Exp(Big2, new(big.Int).Mod(y, k), nil)
}

func tameKangaroo(p, g, b, k *big.Int) (xT, yT *big.Int) {
	N := computeN(k)

	xT = Big0
	yT = new(big.Int).Exp(g, b, p)

	for i := uint64(0); i < N.Uint64(); i++ {
		xT = new(big.Int).Add(xT, f(yT, k))
		yT = new(big.Int).Mod(new(big.Int).Mul(yT, new(big.Int).Exp(g, f(yT, k), p)), p)
	}

	if yT.Cmp(new(big.Int).Exp(g, new(big.Int).Add(b, xT), p)) != 0 {
		fmt.Printf("tame kangaroo did not have expected value")
		return nil, nil
	}

	return
}

// catchKangaroo implements Pollard's kangaroo algorithm.
func catchKangaroo(p, g, y, a, b *big.Int) (m *big.Int, err error) {
	newk := math.Sqrt(float64(p.Uint64()))
	k := new(big.Int).SetInt64(int64(newk))
	// var k = big.NewInt(100000000)
	spew.Dump(k)

	xT, yT := tameKangaroo(p, g, b, k)

	xW := Big0
	yW := y

	// b - a + xT
	bound := new(big.Int).Add(new(big.Int).Sub(b, a), xT)

	for xW.Cmp(bound) < 0 {
		xW = new(big.Int).Add(xW, f(yW, k))
		yW = new(big.Int).Mod(new(big.Int).Mul(yW, new(big.Int).Exp(g, f(yW, k), p)), p)

		if yW.Cmp(yT) == 0 {
			if y.Cmp(new(big.Int).Exp(g, new(big.Int).Add(b, new(big.Int).Sub(xT, xW)), p)) != 0 {
				fmt.Printf("wild  kangaroo did not have expected value")
				return nil, nil
			}
			return new(big.Int).Add(b, new(big.Int).Sub(xT, xW)), nil
		}
	}
	return nil, fmt.Errorf("kangaroo sequences did not intersect")
}

func runDHKangarooAttack(p, g, q, cofactor *big.Int, dh func(*big.Int) []byte, getPublicKey func() *big.Int) (priv *big.Int) {
	var B, R []*big.Int

	factors := findSmallFactors(cofactor)

	for _, order := range factors {
		h, err := findElementOfOrder(order, p)
		if err != nil {
			break
		}
		digest := dh(h)

		b, err := bruteForcedMAC(digest, h, order, p)
		if err != nil {
			fmt.Printf("Error: %v", err)
			break
		}
		B = append(B, b)
		R = append(R, order)
	}

	n, r, err := crt(B, R)
	if err != nil {
		fmt.Printf("Error: %v", err)
		return nil
	}

	h := new(big.Int).Exp(g, r, p)
	y := getPublicKey()

	yN := new(big.Int).Mod(new(big.Int).Mul(y, new(big.Int).Exp(g, new(big.Int).Neg(n), p)), p)

	m, err := catchKangaroo(p, h, yN, Big0, new(big.Int).Div(new(big.Int).Sub(q, Big1), r))
	if err != nil {
		fmt.Printf("Error: %v", err)
		return nil
	}
	return new(big.Int).Sub(n, new(big.Int).Mul(m, r))
}
