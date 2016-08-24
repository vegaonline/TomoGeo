#include "reconstruct.h"
#include "kdnode.h"

template <class T> KDNode<T>::KDNode()
{
    // default constructor
    // ?
}

template <class T> KDNode<T>::KDNode(const T& item, KDNode<T>* nextPtr)
{
    this->data = item;
    this->next = nextPtr;
}

template <class T> void KDNode<T>::InsertAfter(KDNode<T>* p)
{
    p->next = this->next;
    this->next = p;
}

template <class T> KDNode<T>* KDNode<T>::DeleteAfter()
{
    KDNode<T>* tempNode = next;
    if(next != NULL)
        next = next->next;

    return tempNode;
}

template <class T> KDNode<T>* GetNode(const T& item, KDNode<T>* nextPtr = NULL)
{
    KDNode<T>* newNode;
    newNode = new KDNode<T>(item, nextPtr);
    if(newNode == NULL) {
        std::cerr << " Memory allocation failed. " << std::endl;
        exit(1);
    }
    return newNode;
}
